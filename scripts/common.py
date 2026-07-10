"""Shared helpers for the panmap Figure 2 / Figure 3 evaluation.

All panmap commands target the CURRENT build (no --remove-node): leave-one-out is
done by excluding the held-out leaf from the candidate set after scoring, which is
equivalent to pruning the leaf for placement/assembly purposes.
"""
import json
import os
import platform
import re
import subprocess
import sys
import time
from collections import namedtuple
from pathlib import Path


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def _wait_or_kill(proc, upstream, timeout):
    """Wait for a pipeline's terminal consumer `proc`; on timeout kill it and its
    upstream producers so a wedged pipe (e.g. mpileup->ivar) can't hang forever.
    The caller must already have closed its own copy of the inter-process pipes."""
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        for pr in [proc, *upstream]:
            try:
                pr.kill()
            except Exception:
                pass
        proc.wait()


def timed(cmd):
    """Run cmd; return (CompletedProcess, wall_seconds, peak_rss_mb) via /usr/bin/time."""
    flag = "-l" if platform.system() == "Darwin" else "-v"
    t0 = time.monotonic()
    p = subprocess.run(["/usr/bin/time", flag, *cmd], capture_output=True, text=True)
    wall = time.monotonic() - t0
    rss_mb = float("nan")
    for line in p.stderr.splitlines():
        if "maximum resident set size" in line.lower():
            n = float(re.findall(r"\d+", line)[0])
            rss_mb = n / 1e6 if flag == "-l" else n / 1024  # mac: bytes, linux: KB
            break
    return p, wall, rss_mb


# ── panman structure ──────────────────────────────────────────────────────────

_PARENTS_CACHE = {}


def _newick_text(panman, panmanutils):
    """Newick line for a panman, cached on disk (.eval_cache/) so we run the
    (slow) panmanUtils tree export at most once per panman across all jobs."""
    import hashlib
    cache_dir = os.path.join(os.getcwd(), ".eval_cache")
    os.makedirs(cache_dir, exist_ok=True)
    h = hashlib.md5(os.path.abspath(panman).encode()).hexdigest()[:12]
    cache = os.path.join(cache_dir, f"{os.path.basename(panman)}.{h}.nwk")
    if os.path.exists(cache) and os.path.getsize(cache) > 0:
        return open(cache).read()
    nwk = next((l for l in sh([panmanutils, "-I", panman, "-t"]).stdout.splitlines()
                if l.startswith("(")), "")
    if nwk:
        with open(cache, "w") as f:
            f.write(nwk)
    return nwk


def newick_parents(panman, panmanutils):
    """{node_label: parent_label} from the panman's tree (dendropy parse;
    panman internal nodes carry labels, e.g. node_3509). Memoized per panman
    (in-process + on-disk) since it is called once per replicate."""
    key = os.path.abspath(panman)
    if key in _PARENTS_CACHE:
        return _PARENTS_CACHE[key]
    import dendropy
    nwk = _newick_text(panman, panmanutils)
    tree = dendropy.Tree.get(data=nwk, schema="newick", preserve_underscores=True)
    lab = lambda nd: (nd.taxon.label if nd.taxon and nd.taxon.label else nd.label)
    parents = {}
    for nd in tree.preorder_node_iter():
        pl = lab(nd)
        for ch in nd.child_nodes():
            cl = lab(ch)
            if cl is not None:
                parents[cl] = pl
    _PARENTS_CACHE[key] = parents
    return parents


def leaves(parents):
    internal = set(parents.values())
    return [n for n in parents if n not in internal]


def dump_seq(panmap, panman, node, out_fa):
    sh([panmap, panman, "--dump-sequence", node, "-o", out_fa])
    return read_fasta(out_fa)


def ensure_genomes_fa(panmanutils, panman, genomes_fa, samtools):
    """Dump every tip genome to one FASTA once (panmanUtils --fasta) and faidx
    it, so per-node sequence lookups are O(1) instead of reloading the whole
    panman each time (~6.5 s per dump for M. tb). Built lazily under a lock so
    concurrent jobs don't race. panmanUtils writes to <cwd>/info/<pfx>_N.fasta
    (one per PanMAT), which we concatenate into genomes_fa."""
    import fcntl
    import glob
    if os.path.exists(genomes_fa) and os.path.exists(genomes_fa + ".fai"):
        return genomes_fa
    out_dir = os.path.dirname(os.path.abspath(genomes_fa)) or "."
    os.makedirs(out_dir, exist_ok=True)
    with open(genomes_fa + ".buildlock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            if os.path.exists(genomes_fa) and os.path.exists(genomes_fa + ".fai"):
                return genomes_fa
            pu = os.path.abspath(panmanutils) if os.sep in panmanutils else panmanutils
            subprocess.run([pu, "-I", os.path.abspath(panman), "--fasta",
                            "-o", "_allgenomes"], cwd=out_dir,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            parts = sorted(glob.glob(os.path.join(out_dir, "info", "_allgenomes_*.fasta")))
            tmp = genomes_fa + ".tmp"
            with open(tmp, "w") as out:
                for p in parts:
                    with open(p) as f:
                        for line in f:
                            out.write(line)
                    os.remove(p)
            try:
                os.rmdir(os.path.join(out_dir, "info"))
            except OSError:
                pass
            os.replace(tmp, genomes_fa)
            subprocess.run([samtools, "faidx", genomes_fa],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)
    return genomes_fa


def faidx_seq(genomes_fa, node, out_fa, samtools):
    """Extract one record by name from the prebuilt genomes FASTA. Writes out_fa
    and returns the sequence, or "" if the node is absent (e.g. an internal
    node, which is not a tip) so the caller can fall back to a panman dump."""
    p = sh([samtools, "faidx", genomes_fa, node])
    lines = p.stdout.splitlines()
    if p.returncode != 0 or len(lines) < 2 or not lines[0].startswith(">"):
        return ""
    seq = "".join(l for l in lines[1:] if not l.startswith(">"))
    if not seq:
        return ""
    with open(out_fa, "w") as f:
        f.write(p.stdout if p.stdout.endswith("\n") else p.stdout + "\n")
    return seq


def count_reads(fq):
    """Number of reads in a FASTQ (.gz or plain)."""
    import gzip
    op = gzip.open if str(fq).endswith(".gz") else open
    n = 0
    with op(fq, "rt") as f:
        for _ in f:
            n += 1
    return n // 4


def mean_read_len(fq, sample=20000):
    """Mean read length (bp) over up to `sample` reads of a FASTQ (.gz or plain).
    Used to compute real-data coverage from the run's ACTUAL read length instead
    of assuming a fixed 150 bp (real runs range ~44-250 bp)."""
    import gzip
    op = gzip.open if str(fq).endswith(".gz") else open
    tot = n = 0
    with op(fq, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:                      # sequence line
                tot += len(line.rstrip("\n"))
                n += 1
                if n >= sample:
                    break
    return tot / n if n else 0.0


def faidx_get(genomes_fa, node, samtools):
    """Return a node's sequence from the prebuilt genomes FASTA (no file write),
    or '' if absent. Used for ground-truth QC (N-content / length)."""
    p = sh([samtools, "faidx", genomes_fa, node])
    lines = p.stdout.splitlines()
    if p.returncode != 0 or len(lines) < 2 or not lines[0].startswith(">"):
        return ""
    return "".join(l for l in lines[1:] if not l.startswith(">"))


def get_seq(panmanutils, panmap, panman, node, out_fa, genomes_fa, samtools):
    """Sequence for a node: O(1) faidx from the prebuilt tips FASTA, falling
    back to a panman dump for internal/ancestral nodes not present as tips."""
    ensure_genomes_fa(panmanutils, panman, genomes_fa, samtools)
    seq = faidx_seq(genomes_fa, node, out_fa, samtools)
    if seq:
        return seq
    return dump_seq(panmap, panman, node, out_fa)


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


# ── reads: real (ENA download + subsample) and simulated (wgsim) ──────────────

def ena_fastqs(run, out_dir):
    """Download paired FASTQs for an SRA/ENA run. Returns (r1, r2), or
    (None, None) when the run has no paired FASTQ on ENA (some runs marked
    PAIRED only expose a single merged file) or the download fails after
    retries -- the caller then emits a no-call row. Writes via temp+rename and
    a completion marker, so an interrupted download never leaves a partial file
    that a later run would silently trust."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    r1, r2 = out_dir / "raw_1.fastq.gz", out_dir / "raw_2.fastq.gz"
    done = out_dir / "download.ok"
    if done.exists():
        return str(r1), str(r2)
    parts = []
    for _ in range(3):
        rep = sh(["curl", "-sL", "--fail", f"https://www.ebi.ac.uk/ena/portal/api/"
                  f"filereport?accession={run}&result=read_run&fields=fastq_ftp&format=tsv"])
        lines = rep.stdout.strip().splitlines()
        if len(lines) >= 2:
            parts = lines[-1].split("\t")[-1].split(";")
            if parts and parts != [""]:
                break
    u1 = next((u for u in parts if u.endswith("_1.fastq.gz")), None)
    u2 = next((u for u in parts if u.endswith("_2.fastq.gz")), None)
    if not (u1 and u2):
        return None, None                       # no paired data available on ENA
    for url, dst in [(u1, r1), (u2, r2)]:
        tmp = str(dst) + ".part"
        for _ in range(3):
            p = sh(["curl", "-sL", "--fail", "-o", tmp, "https://" + url])
            if p.returncode == 0 and os.path.exists(tmp) and os.path.getsize(tmp) > 0:
                os.replace(tmp, dst)
                break
        else:
            return None, None                    # download failed after retries
    done.touch()
    return str(r1), str(r2)


def ontarget_reads(minimap2, samtools, seqtk, ref_fa, raw1, raw2, cache_dir, threads=8):
    """Keep only read pairs that map to ref_fa (the sample's own GenBank assembly).

    Raw SRA runs are mostly off-target host/background (e.g. RSV can be <1% of a
    library), so subsampling from the whole library makes the coverage axis
    meaningless. Mapping to the sample's own genome first makes downstream
    'coverage' mean depth of that genome. Mapped once per sample, cached and
    locked so concurrent place_real/assemble jobs don't race."""
    import fcntl
    cache_dir = str(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)
    o1 = os.path.join(cache_dir, "ontarget_1.fastq")
    o2 = os.path.join(cache_dir, "ontarget_2.fastq")
    marker = os.path.join(cache_dir, "ontarget.ok")
    with open(os.path.join(cache_dir, "ontarget.lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            if os.path.exists(marker):
                return o1, o2
            names = os.path.join(cache_dir, "ontarget.names")
            # QNAMEs of read pairs with a primary mapping (either mate) to ref_fa.
            mm = subprocess.Popen([minimap2, "-ax", "sr", "-t", str(threads), ref_fa,
                                   raw1, raw2], stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL)
            vw = subprocess.Popen([samtools, "view", "-F", "0x904", "-"],
                                  stdin=mm.stdout, stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL, text=True)
            mm.stdout.close()
            seen = set()
            for line in vw.stdout:
                seen.add(line.split("\t", 1)[0])
            vw.wait(); mm.wait()
            with open(names, "w") as f:
                if seen:
                    f.write("\n".join(sorted(seen)) + "\n")
            # seqtk matches the FASTQ name-token (up to whitespace) == the QNAME,
            # so extracting by name pulls both mates and preserves pairing.
            for src, dst in [(raw1, o1), (raw2, o2)]:
                with open(dst, "w") as f:
                    subprocess.run([seqtk, "subseq", src, names], stdout=f,
                                   stderr=subprocess.DEVNULL)
            open(marker, "w").close()
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)
    return o1, o2


def subsample(seqtk, r1, r2, cov, genome_size, read_len, seed, o1, o2):
    """Deterministic subsample to ~cov x of genome_size (paired). Caps at the
    number of available reads (seqtk returns all when n exceeds the input)."""
    n = max(1, int(cov * genome_size / (2 * read_len)))
    for src, dst in [(r1, o1), (r2, o2)]:
        with open(dst, "w") as f:
            subprocess.run([seqtk, "sample", f"-s{seed}", src, str(n)],
                           stdout=f, stderr=subprocess.DEVNULL)
    return o1, o2


def mutate_genome(ref_fa, mut_rate, indel_frac, seed, out_fa, ts_tv=2.0):
    """Make a novel descendant genome under an HKY85-style substitution model:
    per-site substitution at mut_rate where the new base is drawn ~ pi_j * (kappa if
    the change is a transition else 1). pi is the source genome's empirical base
    composition (so GC-skewed genomes get realistic targets) and kappa is set so the
    expected transition/transversion COUNT ratio equals ts_tv (typ. ~2; ts_tv=None or
    <=0 -> Jukes-Cantor uniform). Indels at mut_rate*indel_frac (50/50 ins/del, 1-9 bp).
    Writes out_fa, returns seq."""
    import random
    from collections import Counter
    rng = random.Random(seed)
    src = read_fasta(ref_fa).upper()
    PUR, PYR = set("AG"), set("CT")
    is_ts = lambda i, j: (i in PUR) == (j in PUR)     # both purine or both pyrimidine
    cnt = Counter(b for b in src if b in "ACGT")
    tot = sum(cnt.values()) or 1
    pi = {b: cnt.get(b, 0) / tot for b in "ACGT"}
    if ts_tv and ts_tv > 0:
        # Calibrate kappa (transition rate multiplier) so the expected ts:tv COUNT
        # ratio == ts_tv, accounting for the per-source-base normalization Z_i (which
        # a closed form ignores). ratio(kappa) is monotone; solve by log-bisection.
        def ratio(kappa):
            num = den = 0.0
            for i in "ACGT":
                Z = sum(pi[j] * (kappa if is_ts(i, j) else 1.0) for j in "ACGT" if j != i)
                if Z <= 0:
                    continue
                for j in "ACGT":
                    if j == i:
                        continue
                    w = pi[i] * pi[j] * (kappa if is_ts(i, j) else 1.0) / Z
                    if is_ts(i, j):
                        num += w
                    else:
                        den += w
            return num / den if den > 0 else float("inf")
        lo, hi = 1e-6, 1e6
        for _ in range(60):
            mid = (lo * hi) ** 0.5
            if ratio(mid) < ts_tv:
                lo = mid
            else:
                hi = mid
        kappa = (lo * hi) ** 0.5
    else:
        kappa = 1.0
    # per-source-base weighted alternative bases (HKY: pi_j * kappa^[transition])
    alts = {i: [(j, pi[j] * (kappa if is_ts(i, j) else 1.0)) for j in "ACGT" if j != i]
            for i in "ACGT"}

    def pick_alt(i):
        ws = alts.get(i)
        if not ws or sum(w for _, w in ws) <= 0:      # non-ACGT / degenerate -> uniform
            return rng.choice([c for c in "ACGT" if c != i])
        r = rng.random() * sum(w for _, w in ws)
        acc = 0.0
        for j, w in ws:
            acc += w
            if r <= acc:
                return j
        return ws[-1][0]

    out = []
    for b in src:
        x = rng.random()
        if x < mut_rate:
            out.append(pick_alt(b))
        elif x < mut_rate + mut_rate * indel_frac:
            if rng.random() < 0.5:                       # insertion
                out.append(b)
                out.append("".join(rng.choice("ACGT") for _ in range(rng.randint(1, 9))))
            # else deletion: emit nothing
        else:
            out.append(b)
    seq = "".join(out)
    with open(out_fa, "w") as f:
        f.write(">descendant\n" + seq + "\n")
    return seq


def sim_reads(wgsim, genome_fa, cov, genome_size, read_len, err, seed, o1, o2):
    """Simulate paired reads from genome_fa at cov (sequencing error only)."""
    n = max(1, int(cov * genome_size / (2 * read_len)))
    sh([wgsim, "-N", str(n), "-1", str(read_len), "-2", str(read_len),
        "-r", "0", "-R", "0", "-e", str(err), "-S", str(seed), genome_fa, o1, o2])
    return o1, o2


# ── distances (parsimony) and assembly accuracy via minimap2 cs ───────────────

# Raw differences of a query vs the interior of a target genome. Substitutions
# and deletions are tallied both as events (indel_events) and as affected bases
# (snps are inherently per-base; del_bases); insertions (query bases absent from
# the target) are tallied separately and never consume target span.
_CsCounts = namedtuple(
    "_CsCounts", "snps indel_events aligned interior del_bases ins_bases ins_events")

# Per-assembly Fig-3 scores. Two accuracy conventions are reported side by side:
#   acc_base  = 100 * matches / interior           -- base-based; this is
#               the "% of the genome correctly reconstructed" the axis label
#               claims, so substitutions, deleted bases and uncovered positions
#               all lower it in proportion to the number of bases affected.
#   acc_event = 100 * cov_frac * (1 - event_err)   -- legacy; each indel counts as
#               a single event regardless of length, so it under-weights large
#               indels (a 100 kb deletion costs the same as a 1 bp one).
AssemblyScores = namedtuple(
    "AssemblyScores",
    "acc_base acc_event snps indel_events del_bases ins_bases aligned interior")
# A method that produced no consensus at all: 0% reconstructed. This is a real
# data point for the coverage-limited baseline, not a dropped no-call.
ZERO_SCORES = AssemblyScores(0.0, 0.0, 0, 0, 0, 0, 0, "")


def _cs_counts(minimap2, target_fa, query_fa, exclude_bp):
    """Base- and event-level differences of QUERY vs the interior of the TARGET
    genome (minimap2 asm20 cs). The first/last exclude_bp of the target are ignored
    for all comparisons -- edge errors are not counted AND those positions are
    removed from the aligned length, so completeness and correctness use the same
    interior. When the query aligns to the target in multiple blocks that overlap on
    target coordinates (e.g. a rearranged or rotated assembly that maps in several
    pieces), each interior target base is counted toward `aligned` only once -- the
    union of covered intervals -- so completeness can never exceed the interior
    length. Returns a _CsCounts (substitutions, indel events, aligned interior
    length, interior length, deleted bases, inserted bases, insertion events)."""
    p = sh([minimap2, "-cx", "asm20", "--cs", target_fa, query_fa])
    snps = indel_events = 0
    del_bases = ins_bases = ins_events = 0
    T = 0
    covered = []                                  # (start, end) target intervals covered by any block
    for line in p.stdout.splitlines():
        col = line.split("\t")
        if len(col) < 12:
            continue
        T = int(col[6])                          # target (truth) genome length
        lo, hi = exclude_bp, T - exclude_bp      # interior in genome coordinates
        cs = next((c[5:] for c in col[12:] if c.startswith("cs:Z:")), None)
        if not cs:
            continue
        pos, i = int(col[7]), 0                   # target start (0-based)
        while i < len(cs):
            c = cs[i]
            if c == ":":                          # identical run
                i += 1; num = ""
                while i < len(cs) and cs[i].isdigit():
                    num += cs[i]; i += 1
                n = int(num)
                a, b = max(pos, lo), min(pos + n, hi)
                if b > a:
                    covered.append((a, b))
                pos += n
            elif c == "*":                        # substitution (1 target base)
                if lo <= pos < hi:
                    covered.append((pos, pos + 1)); snps += 1
                pos += 1; i += 3
            elif c == "+":                        # insertion (query only, no target advance)
                i += 1; ins = 0
                while i < len(cs) and cs[i].isalpha():
                    ins += 1; i += 1
                if lo <= pos < hi:
                    ins_bases += ins; ins_events += 1
            elif c == "-":                        # deletion (d target bases)
                i += 1; d = 0
                while i < len(cs) and cs[i].isalpha():
                    d += 1; i += 1
                a, b = max(pos, lo), min(pos + d, hi)
                if b > a:                         # deleted bases still count as "covered"...
                    covered.append((a, b)); del_bases += b - a   # ...but tracked so acc_base can subtract them
                if lo <= pos < hi:
                    indel_events += 1
                pos += d
            else:
                i += 1
    # Union the covered target intervals so overlapping alignment blocks never
    # double-count: aligned = number of distinct interior target bases covered.
    aligned = 0
    cur_end = -1
    for a, b in sorted(covered):
        if a >= cur_end:
            aligned += b - a; cur_end = b
        elif b > cur_end:
            aligned += b - cur_end; cur_end = b
    return _CsCounts(snps, indel_events, aligned, max(0, T - 2 * exclude_bp),
                     del_bases, ins_bases, ins_events)


def genome_distance(minimap2, a_fa, b_fa, exclude_bp):
    c = _cs_counts(minimap2, a_fa, b_fa, exclude_bp)
    return c.snps + c.indel_events


def genome_distance_raw(minimap2, a_fa, b_fa, exclude_bp):
    """Raw edit distance over the target interior: substitutions + inserted bases +
    deleted bases (indels weighted by LENGTH, not counted as single events)."""
    c = _cs_counts(minimap2, a_fa, b_fa, exclude_bp)
    return c.snps + c.del_bases + c.ins_bases


def assembly_scores(minimap2, consensus_fa, truth_fa, mask_bp):
    """Compare a consensus to the ground-truth genome over the masked interior and
    return both accuracy conventions plus the raw base/event error counts:
      acc_base  = 100 * matches / interior, where matches = aligned - snps -
                  del_bases. This is the "% of the genome correctly
                  reconstructed": substitutions, deleted bases and uncovered
                  positions each lower it by the number of bases they affect.
      acc_event = 100 * cov_frac * (1 - event_err) with event_err counting each
                  indel as one event (legacy Fig-3 metric, kept for comparison).
    Inserted query bases (ins_bases) are reported but excluded from acc_base, which
    measures reconstruction of TRUTH bases only. The mask_bp flanks of the truth
    genome are ignored on both sides so a perfect interior scores 100%."""
    c = _cs_counts(minimap2, truth_fa, consensus_fa, mask_bp)
    if c.aligned == 0 or c.interior == 0:
        return AssemblyScores(0.0, 0.0, c.snps, c.indel_events, c.del_bases,
                              c.ins_bases, c.aligned, c.interior)
    event_err = (c.snps + c.indel_events) / c.aligned
    cov_frac = min(c.aligned / c.interior, 1.0)
    acc_event = 100.0 * cov_frac * (1.0 - event_err)
    matches = c.aligned - c.snps - c.del_bases
    acc_base = 100.0 * max(0, matches) / c.interior
    return AssemblyScores(acc_base, acc_event, c.snps, c.indel_events, c.del_bases,
                          c.ins_bases, c.aligned, c.interior)


def genotyping_accuracy(minimap2, consensus_fa, truth_fa, genome_size, mask_bp):
    """Legacy single-number Fig-3 accuracy (event-based), kept for callers that
    want one value. See assembly_scores for the base-based metric + raw counts.
    genome_size is unused (interior length comes from the truth sequence)."""
    return assembly_scores(minimap2, consensus_fa, truth_fa, mask_bp).acc_event


# ── panmap placement & assembly (current binary) ──────────────────────────────

def build_index(panmap, panman, k, s, l, index, threads):
    p, wall, rss = timed([panmap, panman, "-k", str(k), "-s", str(s), "-l", str(l),
                          "--index-out", index, "--stop", "index", "-f",
                          "--threads", str(threads)])
    return wall, rss


def place(panmap, panman, r1, r2, index, prefix, threads, exclude_self=None,
          force_leaf=True, completeness_weight=0.0):
    """Place reads; return (best_node, all_scores_wall, rss). LOO: drop exclude_self.
    force_leaf restricts placement to observed genomes (leaves); with it off,
    ancestral (internal) nodes are candidates, and completeness_weight>0 re-ranks
    them by logContainment + w*log(len/maxLen) to penalize incomplete nodes."""
    scores = prefix + ".scores.tsv"
    cmd = [panmap, panman, r1, r2, "-i", index, "--stop", "place",
           "--dump-all-scores", scores, "-o", prefix, "--threads", str(threads)]
    if force_leaf:
        cmd.append("--force-leaf")
    if completeness_weight and completeness_weight > 0:
        cmd += ["--completeness-weight", str(completeness_weight)]
    p, wall, rss = timed(cmd)
    rows = []
    if Path(scores).exists():
        with open(scores) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            # panmap places by logContainment (the "[ok] place ... LogC" node) -- score
            # that exact column. Fail loud if it's absent: silently falling back to a
            # different column (e.g. containment) scores a node panmap did NOT place on,
            # which quietly inflates placement/assembly error. Needs a panmap build that
            # emits logContainment in --dump-all-scores (github.com/amkram/panmap PR #80).
            if "logContainment" not in hdr:
                raise RuntimeError(
                    f"{scores}: --dump-all-scores has no 'logContainment' column "
                    f"(header={hdr}); rebuild panmap with the logContainment column "
                    "(github.com/amkram/panmap PR #80) so placement is scored by the "
                    "metric panmap actually uses.")
            mi = hdr.index("logContainment")
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) > mi:
                    rows.append((c[0], float(c[mi])))
    # exclude_self may be a single node id or a collection (its whole duplicate
    # group) so byte-identical twins are also dropped from the candidate set.
    if exclude_self is None:
        excl = set()
    elif isinstance(exclude_self, str):
        excl = {exclude_self}
    else:
        excl = set(exclude_self)
    rows = [r for r in rows if r[0] not in excl]
    best = max(rows, key=lambda r: r[1])[0] if rows else None
    return best, wall, rss


def assemble(panmap, panman, r1, r2, index, ref_node, prefix, threads,
             mutation_spectrum=True):
    """Genotype/assemble against ref_node (panmap's own consensus). Returns
    (consensus_fa, wall, rss). mutation_spectrum=False passes
    --no-mutation-spectrum to disable the substitution-matrix genotype prior."""
    cmd = [panmap, panman, r1, r2, "-i", index, "--reference-node", ref_node,
           "-o", prefix, "--threads", str(threads)]
    if not mutation_spectrum:
        cmd.append("--no-mutation-spectrum")
    p, wall, rss = timed(cmd)
    cons = prefix + ".consensus.fa"
    return (cons if Path(cons).exists() else None), wall, rss


# ── BWA-MEM + iVar baseline (Fig 3, RSV & SARS) ───────────────────────────────

def ensure_bwa_index(bwa, ref_fa):
    """Build the bwa index once under a cross-process lock. A completion
    marker (.bwaok) is written only after `bwa index` finishes, so concurrent jobs
    never race to build the shared reference or read a half-written index."""
    import fcntl
    marker = ref_fa + ".bwaok"
    if Path(marker).exists():
        return
    with open(ref_fa + ".idxlock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            if not Path(marker).exists():
                sh([bwa, "index", ref_fa])
                Path(marker).touch()
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


def impute_to_ref(minimap2, consensus_fa, ref_fa, out_fa):
    """Fill uncovered / N positions of a consensus with the reference base, in
    reference coordinates. iVar emits N wherever depth < -m, which the accuracy
    metric then scores as unreconstructed; imputing to the reference makes the
    single-reference baseline fall back to its reference exactly as panmap falls
    back to its placed reference -- a like-for-like comparison. We start from the
    reference and overlay only the confidently-called (non-N) consensus bases,
    via the minimap2 cs alignment, so indels and long N runs stay coordinate-safe."""
    ref = list(read_fasta(ref_fa))
    p = sh([minimap2, "-cx", "asm20", "--cs", ref_fa, consensus_fa])
    for line in p.stdout.splitlines():
        col = line.split("\t")
        if len(col) < 12:
            continue
        rpos = int(col[7])                       # ref start (0-based)
        cs = next((c[5:] for c in col[12:] if c.startswith("cs:Z:")), None)
        if not cs:
            continue
        i = 0
        while i < len(cs):
            c = cs[i]
            if c == ":":                         # identical run -> keep reference
                i += 1; num = ""
                while i < len(cs) and cs[i].isdigit():
                    num += cs[i]; i += 1
                rpos += int(num)
            elif c == "*":                       # substitution *<refbase><qrybase>
                qb = cs[i + 2]
                if qb.upper() != "N" and rpos < len(ref):
                    ref[rpos] = qb.upper()
                rpos += 1; i += 3
            elif c == "+":                       # insertion in consensus -> ignore (ref coords)
                i += 1
                while i < len(cs) and cs[i].isalpha():
                    i += 1
            elif c == "-":                       # deletion vs consensus -> keep reference
                i += 1
                while i < len(cs) and cs[i].isalpha():
                    rpos += 1; i += 1
            else:
                i += 1
    with open(out_fa, "w") as f:
        f.write(">consensus_imputed\n" + "".join(ref) + "\n")
    return out_fa


def _tree_rss_kb(root):
    """Sum VmRSS (KB) of `root` and all its descendant processes (best-effort)."""
    children = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as f:
                ppid = int(f.read().split()[3])
            children.setdefault(ppid, []).append(int(d))
        except Exception:
            pass
    total, stack = 0, [root]
    while stack:
        pid = stack.pop()
        try:
            with open(f"/proc/{pid}/status") as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        total += int(line.split()[1]); break
        except Exception:
            pass
        stack.extend(children.get(pid, []))
    return total


def timed_tree(cmd, **kw):
    """Run cmd; return (returncode, wall_s, peak_rss_mb) where peak RSS is the max
    summed RSS of the whole process tree, so child tools (gatk/bowtie2/hisat2) count."""
    import threading
    t0 = time.monotonic()
    p = subprocess.Popen(cmd, **kw)
    peak = [0]; stop = threading.Event()

    def sample():
        while not stop.is_set():
            v = _tree_rss_kb(p.pid)
            if v > peak[0]:
                peak[0] = v
            stop.wait(0.15)
    th = threading.Thread(target=sample, daemon=True); th.start()
    p.wait(); stop.set(); th.join(timeout=1)
    return p.returncode, time.monotonic() - t0, peak[0] / 1024.0


def haphpipe_consensus(hp_bin, ref_fa, r1, r2, prefix, bcftools="bcftools"):
    """RSV standard pipeline (Fig 3 Row 1): HaphPIPE finalize_assembly then
    bcftools consensus onto ref_fa. hp_bin is the haphpipe conda env bin dir
    (must supply picard/bowtie2/gatk/samtools). Returns (consensus_fa, wall, rss)."""
    outdir = prefix + "_hp"
    os.makedirs(outdir, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = hp_bin + os.pathsep + env.get("PATH", "")
    _, w_hp, rss_hp = timed_tree([os.path.join(hp_bin, "haphpipe"), "finalize_assembly",
                                  "--ncpu", "1", "--fq1", r1, "--fq2", r2, "--ref_fa", ref_fa,
                                  "--outdir", outdir, "--sample_id", "s"], env=env,
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    vcf = os.path.join(outdir, "final.vcf.gz")
    if not Path(vcf).exists():
        return None, w_hp, rss_hp
    cons = prefix + ".hp.consensus.fa"
    with open(cons, "w") as f:                          # bcftools consensus counted in the timing
        _, w_bc, rss_bc = timed_tree([bcftools, "consensus", "-f", ref_fa, vcf],
                                     stdout=f, stderr=subprocess.DEVNULL)
    return ((cons if Path(cons).exists() and os.path.getsize(cons) else None),
            w_hp + w_bc, max(rss_hp, rss_bc))


def ncbi_sars_consensus(script, ref_fa, hisat2_idx, r1, r2, prefix, custom_filter):
    """SARS-CoV-2 standard pipeline (Fig 3 Row 1): the NCBI SC2VC Illumina
    workflow (Trimmomatic -> HISAT2 -> HaplotypeCaller -> GenotypeGVCFs ->
    VariantFiltration -> LeftAlignAndTrim -> custom_vcf_filter) then bcftools
    consensus, run via scripts/ncbi_sars_pipeline.py. Returns (consensus_fa, wall, rss)."""
    _, wall, rss = timed_tree(["python3", script, ref_fa, hisat2_idx, r1, r2, prefix, custom_filter])
    cons = prefix + ".ncbi.consensus.fa"
    return (cons if Path(cons).exists() and os.path.getsize(cons) else None), wall, rss


def clockwork_consensus(image, ref_dir, ref_fa, r1, r2, prefix, bcftools="bcftools"):
    """M. tuberculosis standard pipeline (Fig 3 Row 1): per-sample Clockwork
    reference_prepare (minimap2 + Cortex reference indexes) then
    variant_call_one_sample (Trimmomatic -> minimap2 -> samtools + Cortex ->
    minos adjudicate), in Docker pinned to one CPU (--cpus 1), then bcftools
    consensus onto ref_fa (H37Rv NC_000962.3). Reference indexing and the
    consensus step are both counted in the wall time; peak container memory is
    polled via `docker stats`. Returns (consensus_fa, wall, peak_mb)."""
    work = os.path.abspath(prefix + "_cw")
    if os.path.exists(work):                     # nuke stale (maybe root-owned) via a throwaway container
        subprocess.run(["docker", "run", "--rm", "-v", f"{os.path.dirname(work)}:/p",
                        "alpine", "rm", "-rf", f"/p/{os.path.basename(work)}"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.makedirs(work)
    reads_dir = os.path.dirname(os.path.abspath(r1))
    ref_src = os.path.dirname(os.path.abspath(ref_fa))
    uid, gid = os.getuid(), os.getgid()
    import threading, hashlib
    tag = hashlib.md5(os.path.abspath(prefix).encode()).hexdigest()[:16]
    base = ["docker", "run", "--rm", "--cpus", "1", "--user", f"{uid}:{gid}", "-w", "/work",
            "-e", "HOME=/work", "-e", "MPLCONFIGDIR=/work/mpl", "-v", f"{work}:/work"]

    def drun(name, extra, cmd):                   # one-CPU docker run, timed with peak-mem polling
        subprocess.run(["docker", "rm", "-f", name], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        t = time.monotonic()
        dp = subprocess.Popen(base + ["--name", name] + extra + [image] + cmd,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pk = [0.0]; stop = threading.Event()

        def mon():
            while not stop.is_set():
                r = subprocess.run(["docker", "stats", "--no-stream", "--format", "{{.MemUsage}}", name],
                                   capture_output=True, text=True)
                mb = _parse_docker_mem(r.stdout)
                if mb > pk[0]:
                    pk[0] = mb
                stop.wait(1.0)
        th = threading.Thread(target=mon, daemon=True); th.start()
        dp.wait(); stop.set(); th.join(timeout=2)
        return time.monotonic() - t, pk[0]

    w_rp, m_rp = drun("cwrp" + tag, ["-v", f"{ref_src}:/rin:ro"],   # reference_prepare, per sample
                      ["clockwork", "reference_prepare", "--outdir", "/work/refprep",
                       f"/rin/{os.path.basename(ref_fa)}"])
    w_vc, m_vc = drun("cw" + tag, ["-v", f"{reads_dir}:/reads:ro"], # variant call vs the fresh reference
                      ["clockwork", "variant_call_one_sample", "--sample_name", "s", "/work/refprep",
                       "/work/out", f"/reads/{os.path.basename(r1)}", f"/reads/{os.path.basename(r2)}"])
    wall = w_rp + w_vc
    rss = max(m_rp, m_vc) or float("nan")
    vcf = os.path.join(work, "out", "final.vcf")
    if not Path(vcf).exists():
        return None, wall, rss
    t = time.monotonic()                         # bcftools consensus counted in the timing
    gz = vcf + ".gz"
    with open(gz, "wb") as f:
        subprocess.run([bcftools, "view", vcf, "-Oz"], stdout=f, stderr=subprocess.DEVNULL)
    subprocess.run([bcftools, "index", "-f", gz], stderr=subprocess.DEVNULL)
    cons = prefix + ".cw.consensus.fa"
    with open(cons, "w") as f:
        subprocess.run([bcftools, "consensus", "-f", ref_fa, gz], stdout=f, stderr=subprocess.DEVNULL)
    wall += time.monotonic() - t
    return (cons if Path(cons).exists() and os.path.getsize(cons) else None), wall, rss


def _parse_docker_mem(s):
    """Parse the first field of `docker stats` MemUsage (e.g. '1.5GiB / 62GiB') to MB."""
    try:
        m = re.match(r"([0-9.]+)\s*([A-Za-z]+)", s.strip().split("/")[0].strip())
        if not m:
            return 0.0
        return float(m.group(1)) * {"gib": 1024, "gb": 1000, "mib": 1, "mb": 1,
                                    "kib": 1 / 1024, "kb": 1 / 1000, "b": 1e-6}.get(m.group(2).lower(), 0)
    except Exception:
        return 0.0


def _clean_ref(ref_fa, out_fa):
    """Write a single-contig copy of ref_fa named '>ref'. iVar consensus calls
    std::stoi on tokens derived from the contig name and crashes
    (std::invalid_argument) when it contains '/', '|' or digits -- which panmap
    leaf labels (e.g. USA/NY.../2022|ON210390.1|2022-03-28) always do. A clean
    name sidesteps the crash; sequence identity (all that the accuracy metric
    reads) is unchanged."""
    seq = []
    with open(ref_fa) as f:
        for line in f:
            if line.startswith(">"):
                if seq:                          # keep only the first record
                    break
            else:
                seq.append(line.strip())
    with open(out_fa, "w") as o:
        o.write(">ref\n" + "".join(seq) + "\n")
    return out_fa


def bwa_ivar(ref_fa, r1, r2, prefix, threads, minimap2=None, bindir="", impute=False):
    """BWA-MEM + iVar consensus. Default leaves iVar's N-gaps at uncovered/low-depth
    positions (the coverage-limited baseline). With impute=True the N positions are
    filled from ref_fa (impute_to_ref). Returns (consensus_fa, wall_seconds)."""
    b = lambda x: os.path.join(bindir, x) if bindir else x
    ref_fa = _clean_ref(ref_fa, prefix + ".cref.fa")   # sanitize contig name for iVar
    bam = prefix + ".bam"
    t0 = time.monotonic()
    sh([b("bwa"), "index", ref_fa])                    # index built per sample, inside timing
    sh([b("samtools"), "faidx", ref_fa])
    _bwaerr = open(prefix + ".bwatime", "w")
    p = subprocess.Popen(["/usr/bin/time", "-v", b("bwa"), "mem", "-t", str(threads), ref_fa, r1, r2],
                         stdout=subprocess.PIPE, stderr=_bwaerr)
    with open(bam, "wb") as bf:
        srt = subprocess.Popen([b("samtools"), "sort", "-o", "-", "-"], stdin=p.stdout,
                               stdout=bf, stderr=subprocess.DEVNULL)
        p.stdout.close()             # parent must drop its read-end (else no EPIPE on crash)
        _wait_or_kill(srt, [p], 3600)
    p.wait()
    # -d 1000: cap pileup depth. Unlimited (-d 0) is pathological on large repeat-rich
    # genomes (e.g. M.tb PE/PPE, IS elements pile to tens of thousands of reads); 1000
    # is far more than consensus needs, so the called sequence is identical.
    mp = subprocess.Popen([b("samtools"), "mpileup", "-aa", "-A", "-d", "1000", "-Q", "20",
                           "--reference", ref_fa, bam], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL)
    iv = subprocess.Popen([b("ivar"), "consensus", "-p", prefix, "-q", "20", "-t", "0.5",
                           "-m", "10"], stdin=mp.stdout, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL)
    mp.stdout.close()                # ivar is the sole reader now: if it exits early,
    _wait_or_kill(iv, [mp], 3600)    # mpileup gets EPIPE instead of hanging
    mp.wait()
    wall = time.monotonic() - t0
    _bwaerr.close()
    rss = float("nan")
    try:
        for _l in open(prefix + ".bwatime"):
            if "Maximum resident set size" in _l:
                rss = float(_l.rsplit(":", 1)[1]) / 1024   # KB -> MB (peak of bwa mem)
    except Exception:
        pass
    raw = prefix + ".fa"
    if not Path(raw).exists():
        return None, wall, rss
    if impute:
        imputed = prefix + ".imputed.fa"
        impute_to_ref(minimap2, raw, ref_fa, imputed)
        return imputed, wall, rss
    return raw, wall, rss

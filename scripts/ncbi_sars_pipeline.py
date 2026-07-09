#!/usr/bin/env python3
"""NCBI SARS-CoV-2 Illumina Variant Calling pipeline (SC2VC, corona_illumina.smk),
run exactly as published, then bcftools consensus. Writes <prefix>.ncbi.consensus.fa.
Usage: ncbi_sars_pipeline.py <ref_fa> <hisat2_idx> <r1> <r2> <prefix> <custom_vcf_filter.py>
Tools are resolved from their known conda envs (Trimmomatic 0.33, HISAT2, GATK4,
Picard, samtools, bcftools, bedtools)."""
import os
import shutil
import subprocess
import sys

ref0, hisat2_idx, r1, r2, prefix, custom_filter = sys.argv[1:7]
SM = "s"                                   # sample name used in the AD filter expression
TRIM = "/home/alex/micromamba/envs/haphpipe/bin/trimmomatic"
HISAT2 = "/home/alex/micromamba/envs/telexpr/bin/hisat2"
HISAT2BUILD = HISAT2 + "-build"
GATK = "/usr/local/bin/gatk"
PICARD = "/home/alex/micromamba/envs/haphpipe/bin/picard"
SAM = "/usr/bin/samtools"
BCF = "/usr/bin/bcftools"
BED = "/usr/bin/bedtools"
DEV = subprocess.DEVNULL


def run(cmd, **kw):
    subprocess.run(cmd, stdout=kw.get("stdout", DEV), stderr=DEV, check=False)


# All reference indexing is rebuilt per sample and counted in this script's runtime:
# copy the reference locally, then faidx + GATK dict + HISAT2 index from scratch (1 thread).
ref = prefix + ".ref.fa"
shutil.copy(ref0, ref)
run([SAM, "faidx", ref])
refdict = os.path.splitext(ref)[0] + ".dict"
run([GATK, "CreateSequenceDictionary", "-R", ref, "-O", refdict])
hisat2_idx = prefix + ".hisat2"
run([HISAT2BUILD, "-p", "1", ref, hisat2_idx])

# 1. Trimmomatic 0.33 (PE)
t1, t2, u1, u2 = (prefix + s for s in (".t1.fq", ".t2.fq", ".u1.fq", ".u2.fq"))
run([TRIM, "PE", "-threads", "1", "-phred33", r1, r2, t1, u1, t2, u2,
     "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"])

# 2. HISAT2 -> -F256 -> sort -> Picard read groups -> index
tmpbam = prefix + ".tmp.bam"
h = subprocess.Popen([HISAT2, "-p", "1", "--no-spliced-alignment", "--no-unal", "-x", hisat2_idx, "-q",
                      "-1", t1, "-2", t2, "-U", f"{u1},{u2}"], stdout=subprocess.PIPE, stderr=DEV)
v = subprocess.Popen([SAM, "view", "-Sb", "-F256", "-"], stdin=h.stdout, stdout=subprocess.PIPE, stderr=DEV)
h.stdout.close()                       # parent drops read-ends so a crashed consumer -> EPIPE, not a hang
with open(tmpbam, "wb") as bf:
    srt = subprocess.Popen([SAM, "sort", "-", "-o", "-"], stdin=v.stdout, stdout=bf, stderr=DEV)
    v.stdout.close()
    srt.wait()
h.wait(); v.wait()
bam = prefix + ".ref.bam"
run([PICARD, "AddOrReplaceReadGroups", f"I={tmpbam}", f"O={bam}",
     "RGID=1", "RGPL=Illumina", "RGPU=NA", f"RGSM={SM}", "RGLB=NA"])
run([SAM, "index", bam])

# 3. HaplotypeCaller (BP_RESOLUTION GVCF)
gvcf = prefix + ".g.vcf.gz"
run([GATK, "HaplotypeCaller", "-R", ref, "-I", bam, "-O", gvcf, "--native-pair-hmm-threads", "1",
     "--minimum-mapping-quality", "10", "--ploidy", "2", "-ERC", "BP_RESOLUTION"])

# 4. GenotypeGVCFs
geno = prefix + ".geno.vcf"
run([GATK, "GenotypeGVCFs", "-R", ref, "-V", gvcf, "-O", geno])

# 5. VariantFiltration (exact 7 filters)
filt = prefix + ".filt.vcf"
run([GATK, "VariantFiltration", "-R", ref, "-V", geno, "-O", filt,
     "--filter-name", "lowAD10", "--filter-expression", f'vc.getGenotype("{SM}").getAD().1 < 10',
     "--filter-name", "lowQUAL100", "--filter-expression", "QUAL < 100",
     "--filter-name", "genomeEnd", "--filter-expression", "POS > 29850",
     "--filter-name", "highFS60", "--filter-expression", "FS >= 60.0",
     "--filter-name", "lowQD2.0", "--filter-expression", "QD < 2.0",
     "--filter-name", "lowReadPosRankSum4.0", "--filter-expression", "ReadPosRankSum < -4.0",
     "--filter-name", "highSOR4.0", "--filter-expression", "SOR >= 4.0"])

# 6. LeftAlignAndTrimVariants
norm = prefix + ".norm.vcf"
run([GATK, "LeftAlignAndTrimVariants", "--verbosity", "ERROR", "--split-multi-allelics",
     "--QUIET", "-R", ref, "-V", filt, "-O", norm])

# 7. bedtools genomecov (per-position depth)
cov = prefix + ".genomecov"
with open(cov, "w") as f:
    run([BED, "genomecov", "-d", "-ibam", bam], stdout=f)

# 8. custom_vcf_filter.py (NCBI's own coverage-based filter)
custom = prefix + ".custom.vcf"
run(["python3", custom_filter, "--c", cov, "--i", norm, "--o", custom])
finalvcf = custom if os.path.exists(custom) else norm

# 9. bcftools consensus on PASS variants
passvcf = prefix + ".pass.vcf.gz"
with open(passvcf, "wb") as f:
    run([BCF, "view", "-f", "PASS", finalvcf, "-Oz"], stdout=f)
run([BCF, "index", passvcf])
cons = prefix + ".ncbi.consensus.fa"
with open(cons, "w") as f:
    run([BCF, "consensus", "-f", ref, passvcf], stdout=f)
print("wrote", cons)

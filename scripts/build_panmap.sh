#!/usr/bin/env bash
# Build panmap + panmanUtils from the latest <ref> of the panmap repo and print the
# panmap binary path (last line of stdout; everything else goes to stderr).
#
#   build_panmap.sh <repo> <build_dir> [ref]
#     repo       path to a panmap git checkout (the eval lives under it; default "..")
#     build_dir  where to keep the build (a detached git worktree, so <repo>'s own
#                branch and working tree are never touched)
#     ref        branch/tag to build (default "main")
#
# Idempotent: fetches <ref>, and rebuilds only when it moved or the binary is missing.
# The eval scores placement by the logContainment column of --dump-all-scores, so a
# build that predates that column (amkram/panmap#80) is rejected here rather than left
# to fail deeper in place().
set -euo pipefail

REPO=$(cd "${1:?repo path}" && pwd)
WT=$(realpath -m "${2:?build dir}")     # absolute: `git -C "$REPO"` would otherwise
REF=${3:-main}                          # resolve a relative worktree path against $REPO
BIN="$WT/build/bin/panmap"

git -C "$REPO" fetch --quiet origin "$REF" 2>/dev/null \
    || echo "build_panmap: fetch failed, using local $REF" >&2
SHA=$(git -C "$REPO" rev-parse "origin/$REF" 2>/dev/null || git -C "$REPO" rev-parse "$REF")

# Point an isolated worktree at SHA (create it on first run).
if [ -e "$WT/.git" ]; then
    git -C "$WT" checkout --quiet --detach "$SHA"
else
    git -C "$REPO" worktree add --quiet --detach "$WT" "$SHA"
fi

if [ "$(cat "$WT/.built_sha" 2>/dev/null || true)" = "$SHA" ] && [ -x "$BIN" ]; then
    echo "$BIN"; exit 0                                  # already current
fi

if ! grep -qF 'containment\tlogContainment' "$WT/src/main.cpp"; then
    echo "build_panmap: $REF has no logContainment column in --dump-all-scores;" >&2
    echo "              merge amkram/panmap#80 into $REF before running the eval." >&2
    exit 1
fi

echo "build_panmap: building panmap @ ${SHA:0:8} (first build compiles deps; be patient)" >&2
cmake -S "$WT" -B "$WT/build" -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=OFF >&2
cmake --build "$WT/build" -j --target panmap panmanUtils >&2
echo "$SHA" > "$WT/.built_sha"
echo "$BIN"

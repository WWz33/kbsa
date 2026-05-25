#!/usr/bin/env bash
# kbsa-anchor-gs2: wrapper that runs GenomeScope2 to estimate kcov,
# then calls `kbsa anchor` with --peak1/--peak2 derived from GS2.
#
# Usage: kbsa-anchor-gs2 <ref.fa> <bulk1_db> <bulk2_db> <bulk1_hist> <bulk2_hist> <out_prefix> [extra kbsa args...]
#
# Required env (or auto-detected from script location):
#   KBSA_BIN          path to kbsa binary  (default: ./build/kbsa)
#   KBSA_GS2_SCRIPT   path to genomescope.R (default: extern/genomescope2/genomescope.R)
#
# Outputs:
#   <hist>.gs2/      GS2 fit directory per bulk (cached; re-used on next call)
#   <hist>.gs2/kcov  one-line cache file with kcov estimate
#   <out_prefix>     standard kbsa anchor outputs

set -euo pipefail

if [[ $# -lt 6 ]]; then
  echo "Usage: $0 <ref.fa> <bulk1_db> <bulk2_db> <bulk1_hist> <bulk2_hist> <out_prefix> [kbsa-args...]" >&2
  exit 1
fi

REF="$1"; BULK1_DB="$2"; BULK2_DB="$3"
BULK1_HIST="$4"; BULK2_HIST="$5"; OUT="$6"
shift 6

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KBSA_BIN="${KBSA_BIN:-${SCRIPT_DIR}/../build/kbsa}"
KBSA_GS2_SCRIPT="${KBSA_GS2_SCRIPT:-${SCRIPT_DIR}/../extern/genomescope2/genomescope.R}"
PLOIDY="${KBSA_PLOIDY:-2}"
KMER="${KBSA_K:-31}"

[[ -x "$KBSA_BIN" ]] || { echo "kbsa binary not executable: $KBSA_BIN" >&2; exit 1; }
[[ -f "$KBSA_GS2_SCRIPT" ]] || { echo "GS2 script not found: $KBSA_GS2_SCRIPT" >&2; exit 1; }

# Run GS2 once per histogram, cache kcov estimate
gs2_kcov() {
  local hist="$1"
  local cache_dir="${hist}.gs2"
  local cache_file="${cache_dir}/kcov"

  if [[ -s "$cache_file" ]]; then
    cat "$cache_file"
    return
  fi

  mkdir -p "$cache_dir"
  Rscript "$KBSA_GS2_SCRIPT" -i "$hist" -o "$cache_dir" -k "$KMER" -p "$PLOIDY" >/dev/null 2>&1 \
    || { echo "GenomeScope2 failed on $hist" >&2; exit 2; }

  local kcov
  kcov=$(awk '/^kmercov/ { print $2; exit }' "$cache_dir/model.txt")
  if [[ -z "$kcov" ]]; then
    echo "Failed to parse kcov from ${cache_dir}/model.txt" >&2
    exit 2
  fi
  printf "%s\n" "$kcov" > "$cache_file"
  printf "%s\n" "$kcov"
}

KCOV1=$(gs2_kcov "$BULK1_HIST")
KCOV2=$(gs2_kcov "$BULK2_HIST")

# diploid (or polyploid) peak = round(ploidy * kcov)
PEAK1=$(awk -v k="$KCOV1" -v p="$PLOIDY" 'BEGIN { printf "%d", k*p + 0.5 }')
PEAK2=$(awk -v k="$KCOV2" -v p="$PLOIDY" 'BEGIN { printf "%d", k*p + 0.5 }')

echo "[kbsa-anchor-gs2] bulk1: kcov=${KCOV1} -> peak1=${PEAK1}" >&2
echo "[kbsa-anchor-gs2] bulk2: kcov=${KCOV2} -> peak2=${PEAK2}" >&2

exec "$KBSA_BIN" anchor \
  --ref "$REF" \
  --bulk1 "$BULK1_DB" --bulk2 "$BULK2_DB" \
  --bulk1-hist "$BULK1_HIST" --bulk2-hist "$BULK2_HIST" \
  --peak1 "$PEAK1" --peak2 "$PEAK2" \
  -o "$OUT" "$@"

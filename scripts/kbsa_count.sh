#!/bin/bash
# kbsa_count.sh — fastp QC + KMC k-mer counting pipeline
# Usage:
#   bash kbsa_count.sh --case-r1 R1.fq.gz --case-r2 R2.fq.gz \
#                      --ctrl-r1 R1.fq.gz --ctrl-r2 R2.fq.gz \
#                      -k 31 -t 20 -o outdir
#   Optional: --parent-r1 / --parent-r2 for trio validation
#   Optional: --skip-qc  to skip fastp (if already cleaned)

set -euo pipefail

# Tool paths (WSL biosoft environment)
FASTP="${FASTP:-/home/ww/biosoft/fastp}"
KMC="${KMC:-/home/ww/miniforge3/bin/kmc}"
KMC_TOOLS="${KMC_TOOLS:-/home/ww/miniforge3/bin/kmc_tools}"

CASE_R1="" CASE_R2="" CTRL_R1="" CTRL_R2=""
PARENT_R1="" PARENT_R2=""
K=31 THREADS=20 MEM=16 OUTDIR=""
SKIP_QC=0 CI=2

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-r1) CASE_R1="$2"; shift 2;;
    --case-r2) CASE_R2="$2"; shift 2;;
    --ctrl-r1) CTRL_R1="$2"; shift 2;;
    --ctrl-r2) CTRL_R2="$2"; shift 2;;
    --parent-r1) PARENT_R1="$2"; shift 2;;
    --parent-r2) PARENT_R2="$2"; shift 2;;
    -k) K="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -m) MEM="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    --ci) CI="$2"; shift 2;;
    --skip-qc) SKIP_QC=1; shift;;
    *) echo "Unknown: $1"; exit 1;;
  esac
done

if [ -z "$CASE_R1" ] || [ -z "$CTRL_R1" ] || [ -z "$OUTDIR" ]; then
  echo "Usage: $0 --case-r1 R1.fq --case-r2 R2.fq --ctrl-r1 R1.fq --ctrl-r2 R2.fq -o outdir"
  exit 1
fi

# Validate tools
for tool in "$FASTP" "$KMC" "$KMC_TOOLS"; do
  if [ ! -x "$tool" ]; then
    echo "ERROR: $tool not found or not executable"; exit 1
  fi
done

TMPDIR="$OUTDIR/tmp"
QCDIR="$OUTDIR/qc"
mkdir -p "$OUTDIR" "$TMPDIR" "$QCDIR"

echo "================================================================"
echo "[kbsa count] k=$K, ci=$CI, threads=$THREADS, mem=${MEM}GB"
echo "[kbsa count] Output: $OUTDIR"
echo "================================================================"

# --- Phase 1: fastp QC ---
run_fastp() {
  local prefix="$1" r1="$2" r2="$3"
  local out1="$QCDIR/${prefix}_clean_R1.fq.gz"
  local out2="$QCDIR/${prefix}_clean_R2.fq.gz"

  if [ "$SKIP_QC" -eq 1 ]; then
    echo "[fastp] Skipped for $prefix (--skip-qc)"
    echo "$r1" > "$TMPDIR/${prefix}.lst"
    [ -n "$r2" ] && echo "$r2" >> "$TMPDIR/${prefix}.lst"
    return
  fi

  echo "[fastp] Processing $prefix..."
  "$FASTP" -i "$r1" -I "$r2" \
    -o "$out1" -O "$out2" \
    --dedup \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread $(( THREADS > 8 ? 8 : THREADS )) \
    --json "$QCDIR/${prefix}_fastp.json" \
    --html "$QCDIR/${prefix}_fastp.html" \
    > "$QCDIR/${prefix}_fastp.log" 2>&1

  echo "$out1" > "$TMPDIR/${prefix}.lst"
  echo "$out2" >> "$TMPDIR/${prefix}.lst"
  echo "[fastp] $prefix done. See $QCDIR/${prefix}_fastp.json"
}

run_fastp "case" "$CASE_R1" "$CASE_R2" &
PID_QC1=$!
run_fastp "ctrl" "$CTRL_R1" "$CTRL_R2" &
PID_QC2=$!

if [ -n "$PARENT_R1" ]; then
  run_fastp "parent" "$PARENT_R1" "$PARENT_R2" &
  PID_QC3=$!
fi

wait $PID_QC1 || { echo "ERROR: fastp case failed"; exit 1; }
wait $PID_QC2 || { echo "ERROR: fastp ctrl failed"; exit 1; }
[ -n "${PID_QC3:-}" ] && { wait $PID_QC3 || { echo "ERROR: fastp parent failed"; exit 1; }; }

echo "[fastp] All QC complete."
echo ""

# --- Phase 2: KMC counting (parallel) ---
run_kmc() {
  local prefix="$1"
  local kmc_tmp="$TMPDIR/${prefix}_kmc_tmp"
  mkdir -p "$kmc_tmp"
  echo "[kmc] Counting $prefix pool..."
  "$KMC" -k$K -t$THREADS -m$MEM -ci$CI -cs65535 \
    -fq @"$TMPDIR/${prefix}.lst" \
    "$OUTDIR/${prefix}_db" "$kmc_tmp" \
    > "$TMPDIR/${prefix}_kmc.log" 2>&1
  rm -rf "$kmc_tmp"
  echo "[kmc] $prefix done: $(grep 'No. of unique counted' "$TMPDIR/${prefix}_kmc.log" 2>/dev/null || echo 'see log')"
}

run_kmc "case" || { echo "ERROR: kmc case failed"; cat "$TMPDIR/case_kmc.log"; exit 1; }
run_kmc "ctrl" || { echo "ERROR: kmc ctrl failed"; cat "$TMPDIR/ctrl_kmc.log"; exit 1; }

if [ -n "$PARENT_R1" ]; then
  run_kmc "parent" || { echo "ERROR: kmc parent failed"; exit 1; }
fi

echo "[kmc] All counting complete."
echo ""

# --- Phase 3: Histograms ---
echo "[hist] Generating histograms..."
"$KMC_TOOLS" transform "$OUTDIR/case_db" histogram "$OUTDIR/case_hist.txt"
"$KMC_TOOLS" transform "$OUTDIR/ctrl_db" histogram "$OUTDIR/ctrl_hist.txt"
[ -f "$OUTDIR/parent_db.kmc_pre" ] && \
  "$KMC_TOOLS" transform "$OUTDIR/parent_db" histogram "$OUTDIR/parent_hist.txt"

# --- Summary ---
echo ""
echo "================================================================"
echo "[kbsa count] Pipeline complete!"
echo ""
echo "  Case DB:  $OUTDIR/case_db  ($(wc -l < "$OUTDIR/case_hist.txt") histogram bins)"
echo "  Ctrl DB:  $OUTDIR/ctrl_db  ($(wc -l < "$OUTDIR/ctrl_hist.txt") histogram bins)"
[ -f "$OUTDIR/parent_hist.txt" ] && \
  echo "  Parent DB: $OUTDIR/parent_db"
echo ""
echo "  QC reports: $QCDIR/"
echo ""
echo "[Next] Run kbsa scoring:"
echo "  kbsa score --case $OUTDIR/case_db --ctrl $OUTDIR/ctrl_db \\"
echo "             --case-hist $OUTDIR/case_hist.txt --ctrl-hist $OUTDIR/ctrl_hist.txt \\"
echo "             -o $OUTDIR/kbsa_scored.tsv"
echo "================================================================"

#!/bin/bash
set -e

KMC=/home/ww/miniforge3/bin/kmc
KMC_TOOLS=/home/ww/miniforge3/bin/kmc_tools
KBSA=/home/ww/biosoft/kbsa
INDIR=/home/ww/kqtl_data/kqtl_cucumber
OUTDIR=/home/ww/kqtl_data/results/cucumber_csdw3

# Pool swap: kbsa Case = WT = P2, kbsa Ctrl = Mutant = P1
CASE_R1=$INDIR/P2_clean_r1.fastq.gz
CASE_R2=$INDIR/P2_clean_r2.fastq.gz
CTRL_R1=$INDIR/P1_clean_r1.fastq.gz
CTRL_R2=$INDIR/P1_clean_r2.fastq.gz

mkdir -p "$OUTDIR"

echo '[cucumber] Pool assignment: Case=P2(WT), Ctrl=P1(Mutant)'

echo '[cucumber] Counting case pool (P2/WT)...'
mkdir -p /dev/shm/kmc_case_tmp
printf '%s\n' "$CASE_R1" "$CASE_R2" > /tmp/cuc_case.lst
$KMC -k31 -ci2 -cs65535 -t10 -m8 -fq @/tmp/cuc_case.lst "$OUTDIR/case_db" /dev/shm/kmc_case_tmp
rm -rf /dev/shm/kmc_case_tmp

echo '[cucumber] Counting ctrl pool (P1/Mutant)...'
mkdir -p /dev/shm/kmc_ctrl_tmp
printf '%s\n' "$CTRL_R1" "$CTRL_R2" > /tmp/cuc_ctrl.lst
$KMC -k31 -ci2 -cs65535 -t10 -m8 -fq @/tmp/cuc_ctrl.lst "$OUTDIR/ctrl_db" /dev/shm/kmc_ctrl_tmp
rm -rf /dev/shm/kmc_ctrl_tmp

echo '[cucumber] Generating histograms...'
$KMC_TOOLS transform "$OUTDIR/case_db" histogram "$OUTDIR/case_hist.txt"
$KMC_TOOLS transform "$OUTDIR/ctrl_db" histogram "$OUTDIR/ctrl_hist.txt"

echo '[cucumber] Running kbsa score...'
$KBSA score \
  --case "$OUTDIR/case_db" \
  --ctrl "$OUTDIR/ctrl_db" \
  --case-hist "$OUTDIR/case_hist.txt" \
  --ctrl-hist "$OUTDIR/ctrl_hist.txt" \
  -o "$OUTDIR/kbsa_scored.tsv"

echo '[cucumber] ALL DONE'
wc -l "$OUTDIR/kbsa_scored.tsv"

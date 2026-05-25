#!/bin/bash
set -euo pipefail

# =============================================================================
# BSA-seq Simulation Pipeline
# Generates synthetic 30x F2 BSA-seq data and tests kbsa detection
#
# Requirements: python3 (numpy), art_illumina, kmc, samtools, seqkit, kbsa
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
KBSA_DIR="$(dirname "$SCRIPT_DIR")"
KBSA_BIN="${KBSA_DIR}/build/kbsa"
OUTDIR="${SCRIPT_DIR}/output"
THREADS=4

# Simulation parameters
CHROM_LENGTH=26000000    # 26Mb (cucumber chr1 scale)
SNP_DENSITY=0.002        # 1 SNP per 500bp between parents
CAUSAL_POS=13000000      # Middle of chromosome
N_F2=2000
BULK_SIZE=40
COVERAGE=30
SEED=42

# kbsa parameters
KMER=31
PEAK_TAU=2.0

echo "============================================"
echo " BSA-seq Simulation Pipeline"
echo "============================================"
echo "Chrom: ${CHROM_LENGTH}bp, SNP density: ${SNP_DENSITY}"
echo "Causal: pos ${CAUSAL_POS}"
echo "F2: ${N_F2} individuals, bulk size: ${BULK_SIZE}"
echo "Coverage: ${COVERAGE}x PE150"
echo "============================================"

# Step 1: Generate F2 bulks
echo ""
echo "=== STEP 1: Generate F2 population and select bulks ==="
python3 "${SCRIPT_DIR}/make_f2_bulks.py" \
    --chrom-length ${CHROM_LENGTH} \
    --snp-density ${SNP_DENSITY} \
    --causal-pos ${CAUSAL_POS} \
    --n-f2 ${N_F2} \
    --bulk-size ${BULK_SIZE} \
    --seed ${SEED} \
    --outdir "${OUTDIR}"

# Step 2: Generate reads with ART
echo ""
echo "=== STEP 2: Generate Illumina PE150 reads (${COVERAGE}x) ==="

GENOME_SIZE=${CHROM_LENGTH}
HAPS=$((BULK_SIZE * 2))  # diploid = 2 haplotypes per individual
# coverage per haplotype = total_coverage / n_haplotypes
# ART uses fold coverage directly
HAP_COV=$(python3 -c "print(f'{${COVERAGE}/${HAPS}:.4f}')")

echo "Per-haplotype coverage: ${HAP_COV}x (${HAPS} haplotypes total)"

mkdir -p "${OUTDIR}/reads_tmp/bulk1" "${OUTDIR}/reads_tmp/bulk2"

generate_reads_for_bulk() {
    local bulk_name=$1
    local seed_offset=$2
    local hap_dir="${OUTDIR}/haps/${bulk_name}"
    local reads_dir="${OUTDIR}/reads_tmp/${bulk_name}"
    local count=0

    for fa in "${hap_dir}"/*.fa; do
        base=$(basename "$fa" .fa)
        art_illumina \
            -ss HSXt \
            -i "$fa" \
            -p \
            -l 150 \
            -f "${HAP_COV}" \
            -m 350 \
            -s 50 \
            -rs $((SEED + seed_offset + count)) \
            -o "${reads_dir}/${base}_" \
            -q &>/dev/null
        count=$((count + 1))
    done

    echo "  ${bulk_name}: generated reads from ${count} haplotypes"
}

echo "Generating bulk1 reads..."
generate_reads_for_bulk "bulk1" 100

echo "Generating bulk2 reads..."
generate_reads_for_bulk "bulk2" 200

# Step 3: Concatenate reads per bulk
echo ""
echo "=== STEP 3: Merge reads ==="

cat "${OUTDIR}/reads_tmp/bulk1/"*_1.fq > "${OUTDIR}/bulk1_R1.fq"
cat "${OUTDIR}/reads_tmp/bulk1/"*_2.fq > "${OUTDIR}/bulk1_R2.fq"
cat "${OUTDIR}/reads_tmp/bulk2/"*_1.fq > "${OUTDIR}/bulk2_R1.fq"
cat "${OUTDIR}/reads_tmp/bulk2/"*_2.fq > "${OUTDIR}/bulk2_R2.fq"

echo "  bulk1: $(wc -l < "${OUTDIR}/bulk1_R1.fq" | awk '{print $1/4}') read pairs"
echo "  bulk2: $(wc -l < "${OUTDIR}/bulk2_R1.fq" | awk '{print $1/4}') read pairs"

# Step 4: Add PCR duplicates (~10%)
echo ""
echo "=== STEP 4: Add PCR duplicates (10%) ==="

for bulk in bulk1 bulk2; do
    for read in R1 R2; do
        fq="${OUTDIR}/${bulk}_${read}.fq"
        seqkit sample -p 0.10 -s $((SEED + 300)) "$fq" > "${fq}.dup"
        cat "$fq" "${fq}.dup" > "${fq}.withdup"
        mv "${fq}.withdup" "$fq"
        rm "${fq}.dup"
    done
done

echo "  Final bulk1: $(wc -l < "${OUTDIR}/bulk1_R1.fq" | awk '{print $1/4}') read pairs"
echo "  Final bulk2: $(wc -l < "${OUTDIR}/bulk2_R1.fq" | awk '{print $1/4}') read pairs"

# Step 5: Build KMC databases
echo ""
echo "=== STEP 5: Build KMC k-mer databases ==="

mkdir -p "${OUTDIR}/kmc_tmp"

# KMC from paired FASTQ
for bulk in bulk1 bulk2; do
    # Create file list for KMC
    echo "${OUTDIR}/${bulk}_R1.fq" > "${OUTDIR}/${bulk}_files.txt"
    echo "${OUTDIR}/${bulk}_R2.fq" >> "${OUTDIR}/${bulk}_files.txt"

    kmc -k${KMER} -t${THREADS} -ci2 -cs65535 \
        @"${OUTDIR}/${bulk}_files.txt" \
        "${OUTDIR}/${bulk}_kmc" \
        "${OUTDIR}/kmc_tmp" \
        2>&1 | grep -E "^(Total|Unique)"

    # Sort KMC database (required for kbsa)
    kmc_tools transform "${OUTDIR}/${bulk}_kmc" sort "${OUTDIR}/${bulk}_sorted" \
        2>&1 | tail -1

    # Histogram (required for kbsa normalization)
    kmc_tools transform "${OUTDIR}/${bulk}_kmc" histogram "${OUTDIR}/${bulk}_hist.txt" \
        2>&1 | tail -1
done

# Step 6: Index reference for kbsa
echo ""
echo "=== STEP 6: Index reference ==="

samtools faidx "${OUTDIR}/ref.fa"

# Step 7: Run kbsa
echo ""
echo "=== STEP 7: Run kbsa anchor (excess-mass tau=${PEAK_TAU}) ==="

"${KBSA_BIN}" anchor \
    --ref "${OUTDIR}/ref.fa" \
    --bulk1 "${OUTDIR}/bulk1_sorted" \
    --bulk2 "${OUTDIR}/bulk2_sorted" \
    --bulk1-hist "${OUTDIR}/bulk1_hist.txt" \
    --bulk2-hist "${OUTDIR}/bulk2_hist.txt" \
    --peak \
    --peak-tau ${PEAK_TAU} \
    --threads ${THREADS} \
    -o "${OUTDIR}/sim_result.bed" \
    2>&1

# Step 8: Evaluate results
echo ""
echo "=== STEP 8: Evaluate detection ==="
echo ""
echo "Causal SNP position: ${CAUSAL_POS} (chr1:${CAUSAL_POS})"
echo "Expected 1Mb window: $((CAUSAL_POS / 1000000))"
echo ""
echo "Top 10 windows by excess-mass score:"
sort -k5 -rn "${OUTDIR}/sim_result.bed" | head -10
echo ""

# Check if causal window is in top N
CAUSAL_MB=$((CAUSAL_POS / 1000000))
CAUSAL_START=$((CAUSAL_MB * 1000000))
CAUSAL_END=$(((CAUSAL_MB + 1) * 1000000))

RANK=$(sort -k5 -rn "${OUTDIR}/sim_result.bed" | \
    awk -v s=${CAUSAL_START} -v e=${CAUSAL_END} \
    'BEGIN{rank=0} {rank++; if($2>=s && $2<e){print rank; found=1; exit}} END{if(!found) print "NOT_FOUND"}')

echo "Causal window rank: #${RANK}"
echo ""

if [ "$RANK" -le 3 ]; then
    echo "PASS: Causal SNP detected in top 3"
elif [ "$RANK" -le 10 ]; then
    echo "MARGINAL: Causal SNP in top 10 (rank #${RANK})"
else
    echo "FAIL: Causal SNP not in top 10 (rank #${RANK})"
fi

# Cleanup temp
rm -rf "${OUTDIR}/reads_tmp" "${OUTDIR}/kmc_tmp"
rm -f "${OUTDIR}"/*_files.txt

echo ""
echo "=== Pipeline complete ==="
echo "Results: ${OUTDIR}/sim_result.bed"
echo "Truth:   ${OUTDIR}/variants_truth.tsv"

#!/bin/bash
set -e

# Brapa ground truth validation
# Target: Bra032670, A09:38876601-38881010, PAV (gene absent in mutant)
# Expect: CONTROL-enriched (WT-enriched) k-mers cluster at this locus

SCORED=/home/ww/kqtl_data/results/brapa/kbsa_scored.tsv
REF=/mnt/f/kqtl/data/brapa_true/Brapa_sequence_v1.5.fa
OUTDIR=/home/ww/kqtl_data/results/brapa/validation
TARGET_CHR="A09"
TARGET_START=38876601
TARGET_END=38881010

mkdir -p "$OUTDIR"

echo "=== Brapa Ground Truth Validation ==="
echo "Target: $TARGET_CHR:$TARGET_START-$TARGET_END (Bra032670 PAV)"
echo ""

# Extract top CONTROL-enriched k-mers (g_score >= 50)
echo "[1] Extracting top CONTROL-enriched k-mers (g>=50)..."
awk -F'\t' '$2=="CONTROL" && $4>=50 {print ">"NR"\n"$1}' "$SCORED" > "$OUTDIR/ctrl_top_kmers.fa"
CTRL_COUNT=$(grep -c '^>' "$OUTDIR/ctrl_top_kmers.fa")
echo "    CONTROL-enriched k-mers (g>=50): $CTRL_COUNT"

# Extract top CASE-enriched k-mers (g>=50) for comparison
awk -F'\t' '$2=="CASE" && $4>=50 {print ">"NR"\n"$1}' "$SCORED" > "$OUTDIR/case_top_kmers.fa"
CASE_COUNT=$(grep -c '^>' "$OUTDIR/case_top_kmers.fa")
echo "    CASE-enriched k-mers (g>=50): $CASE_COUNT"

# Map with bowtie2 (if available) or bwa
if command -v bowtie2 &>/dev/null; then
    MAPPER="bowtie2"
elif command -v bwa &>/dev/null; then
    MAPPER="bwa"
else
    echo "ERROR: No mapper (bowtie2/bwa) found"
    exit 1
fi
echo "[2] Using mapper: $MAPPER"

# Use bwa mem for short sequences
echo "[3] Mapping CONTROL k-mers to reference..."
bwa mem -t 10 -k 31 -T 30 "$REF" "$OUTDIR/ctrl_top_kmers.fa" 2>/dev/null | \
    awk -v chr="$TARGET_CHR" -v s="$TARGET_START" -v e="$TARGET_END" '
    BEGIN {total=0; on_target=0; on_chr=0}
    /^@/ {next}
    {
        total++
        if ($3 == chr) {
            on_chr++
            pos = $4
            if (pos >= s && pos <= e) on_target++
        }
    }
    END {
        printf "    Total mapped: %d\n", total
        printf "    On %s: %d (%.2f%%)\n", chr, on_chr, 100.0*on_chr/total
        printf "    On target region: %d (%.4f%%)\n", on_target, 100.0*on_target/total
        printf "    Enrichment vs random: %.1fx\n", (on_target*1.0/total) / ((e-s)*1.0/400000000)
    }' > "$OUTDIR/ctrl_mapping_stats.txt"
cat "$OUTDIR/ctrl_mapping_stats.txt"

echo ""
echo "[4] Mapping CASE k-mers to reference (negative control)..."
bwa mem -t 10 -k 31 -T 30 "$REF" "$OUTDIR/case_top_kmers.fa" 2>/dev/null | \
    awk -v chr="$TARGET_CHR" -v s="$TARGET_START" -v e="$TARGET_END" '
    BEGIN {total=0; on_target=0; on_chr=0}
    /^@/ {next}
    {
        total++
        if ($3 == chr) {
            on_chr++
            pos = $4
            if (pos >= s && pos <= e) on_target++
        }
    }
    END {
        printf "    Total mapped: %d\n", total
        printf "    On %s: %d (%.2f%%)\n", chr, on_chr, 100.0*on_chr/total
        printf "    On target region: %d (%.4f%%)\n", on_target, 100.0*on_target/total
        printf "    Enrichment vs random: %.1fx\n", (on_target*1.0/total) / ((e-s)*1.0/400000000)
    }' > "$OUTDIR/case_mapping_stats.txt"
cat "$OUTDIR/case_mapping_stats.txt"

echo ""
echo "=== DONE ==="

#!/bin/bash
FILE=/home/ww/kqtl_data/results/brapa/kbsa_scored.tsv

echo "=== G-score distribution ==="
awk -F'\t' 'NR>1 {if($4>=100) c100++; else if($4>=50) c50++; else c20++} END{printf "g>=100: %d\ng>=50: %d\ng<50: %d\n", c100, c50, c20}' "$FILE"

echo ""
echo "=== Top 10 CASE-enriched by kai_reg ==="
awk -F'\t' '$2=="CASE"' "$FILE" | sort -t$'\t' -k3 -nr | head -10

echo ""
echo "=== Top 10 CONTROL-enriched by kai_reg ==="
awk -F'\t' '$2=="CONTROL"' "$FILE" | sort -t$'\t' -k3 -nr | head -10

echo ""
echo "=== kai_reg distribution ==="
awk -F'\t' 'NR>1 {if($3>=0.9) c9++; else if($3>=0.8) c8++; else c0++} END{printf "kai>=0.9: %d\nkai>=0.8: %d\nkai<0.8: %d\n", c9, c8, c0}' "$FILE"

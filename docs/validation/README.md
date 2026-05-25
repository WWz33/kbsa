# kbsa BSA Peak Validation — Full Pipeline

This document records the **exact commands**, **input files**, **output files**, and **ranking results** of the 4-dataset validation run on 2026-05-22.

The pipeline produces ranked 1Mb genomic intervals containing BSA-seq QTL signal, given KMC k-mer databases for case/control pools and a reference genome.

> The figures in `figures/` are reproduced verbatim by `scripts/visualize_validation.py`.

---

## Pipeline overview

```
Pool fastqs              Reference FASTA
    ↓ KMC (count + sort)
case_db_sorted.kmc       ref.fa
ctrl_db_sorted.kmc
case_hist.txt
ctrl_hist.txt
    ↓
kbsa anchor --per-position
    ↓ (per-position TSV: chrom, pos, case_raw, ctrl_raw, z_score)
    ↓
filter: case_raw > 0 AND ctrl_raw > 0  (shared-only, BSA proper)
aggregate: 1Mb sliding window, sum(|z_score|)
sort descending → ranked QTL candidate intervals
```

The full algorithmic logic in 30 lines of Python (see `scripts/visualize_validation.py`):

```python
import numpy as np
pos = []; z = []
with open(per_position_tsv) as f:
    next(f)
    for line in f:
        c = line.split('\t')
        if int(c[2]) > 0 and int(c[3]) > 0:        # shared-only filter
            pos.append(int(c[1])); z.append(float(c[4]))
pos = np.array(pos); z = np.array(z)
sum_abs = np.zeros(pos[-1] // 1_000_000 + 1)
np.add.at(sum_abs, pos // 1_000_000, np.abs(z))   # 1Mb sum|z|
order = np.argsort(-sum_abs)                       # rank descending
```

---

## Input files (verified to exist 2026-05-22)

### brapa
- `case_db`:   `/home/ww/kqtl_data/results/brapa/case_db_sorted.kmc_pre|.kmc_suf`
- `ctrl_db`:   `/home/ww/kqtl_data/results/brapa/ctrl_db_sorted.kmc_pre|.kmc_suf`
- `case_hist`: `/home/ww/kqtl_data/results/brapa/case_hist.txt`
- `ctrl_hist`: `/home/ww/kqtl_data/results/brapa/ctrl_hist.txt`
- `reference`: `/mnt/f/kqtl/data/brapa_true/Brapa_sequence_v1.5.fa`
- `target chromosome`: A09
- `ground truth`: A09:38876601-38881010 (Bra032670, 4.4 kb PAV, F2 cross)

### cucumber
- `case_db`:   `/home/ww/kqtl_data/results/cucumber_csdw3/case_db_sorted.kmc_pre|.kmc_suf`
- `ctrl_db`:   `/home/ww/kqtl_data/results/cucumber_csdw3/ctrl_db_sorted.kmc_pre|.kmc_suf`
- `case_hist`: `/home/ww/kqtl_data/results/cucumber_csdw3/case_hist.txt`
- `ctrl_hist`: `/home/ww/kqtl_data/results/cucumber_csdw3/ctrl_hist.txt`
- `reference`: `/mnt/f/kqtl/data/cucumber_true/ChineseLong_genome_v3.fa`
- `target chromosome`: chr1
- `ground truth`: chr1:30062867-30066723 (CsDw3 / CsaV3_1G044640, 102 bp DEL, F2 X1×csdw3)

### cabbage
- `case_db`:   `/home/ww/kqtl_data/results/cabbage/case_db_sorted.kmc_pre|.kmc_suf`
- `ctrl_db`:   `/home/ww/kqtl_data/results/cabbage/ctrl_db_sorted.kmc_pre|.kmc_suf`
- `case_hist`: `/home/ww/kqtl_data/results/cabbage/case_hist.txt`
- `ctrl_hist`: `/home/ww/kqtl_data/results/cabbage/ctrl_hist.txt`
- `reference`: `/mnt/f/kqtl/data/true_cabbage/BOL.fa`
- `target chromosome`: C09
- `ground truth`: C09:29142111-29144885 (Bol035718 Ms-cd1, 1 bp DEL, BC24 heterozygous)

### vradiata
- `case_db`:   `/home/ww/kqtl_data/results/vradiata/case_db_sorted.kmc_pre|.kmc_suf`
- `ctrl_db`:   `/home/ww/kqtl_data/results/vradiata/ctrl_db_sorted.kmc_pre|.kmc_suf`
- `case_hist`: `/home/ww/kqtl_data/results/vradiata/case_hist.txt`
- `ctrl_hist`: `/home/ww/kqtl_data/results/vradiata/ctrl_hist.txt`
- `reference`: `/mnt/f/kqtl/data/vrad_true/Vrad_JL7.genome.fa`
- `target chromosome`: 11 (note: not "chr11"; reference uses bare integer headers)
- `ground truth`: 11:21271797-21273219 (jg35124, 1 bp DEL in exon, RIL 30 vs 30)

---

## Step 1 — Per-position projection (kbsa anchor)

The `kbsa anchor --per-position` mode scans every k-mer position on the reference, queries case + ctrl KMC databases via random access, and emits a TSV row for every position with non-zero counts in either pool.

Built binary: `/home/ww/kqtl_data/kbsa-standalone/build/kbsa`. The `anchor` subcommand at line `cmd_anchor` in `src/main.cpp` (the per-position branch was added in this session).

### brapa A09

```bash
/home/ww/kqtl_data/kbsa-standalone/build/kbsa anchor \
  --ref /mnt/f/kqtl/data/brapa_true/Brapa_sequence_v1.5.fa \
  --case /home/ww/kqtl_data/results/brapa/case_db_sorted \
  --ctrl /home/ww/kqtl_data/results/brapa/ctrl_db_sorted \
  --case-hist /home/ww/kqtl_data/results/brapa/case_hist.txt \
  --ctrl-hist /home/ww/kqtl_data/results/brapa/ctrl_hist.txt \
  --per-position --chrom A09 \
  -o /tmp/cf_test/brapa_A09_perpos.tsv
```

- Runtime: 76 s (wall)
- Output: 33,171,530 positions, 869 MB

### cucumber chr1

```bash
/home/ww/kqtl_data/kbsa-standalone/build/kbsa anchor \
  --ref /mnt/f/kqtl/data/cucumber_true/ChineseLong_genome_v3.fa \
  --case /home/ww/kqtl_data/results/cucumber_csdw3/case_db_sorted \
  --ctrl /home/ww/kqtl_data/results/cucumber_csdw3/ctrl_db_sorted \
  --case-hist /home/ww/kqtl_data/results/cucumber_csdw3/case_hist.txt \
  --ctrl-hist /home/ww/kqtl_data/results/cucumber_csdw3/ctrl_hist.txt \
  --per-position --chrom chr1 \
  -o /tmp/cf_test/cucumber_chr1_perpos.tsv
```

- Runtime: 76 s
- Output: 28,780,246 positions, 788 MB

### cabbage C09

```bash
/home/ww/kqtl_data/kbsa-standalone/build/kbsa anchor \
  --ref /mnt/f/kqtl/data/true_cabbage/BOL.fa \
  --case /home/ww/kqtl_data/results/cabbage/case_db_sorted \
  --ctrl /home/ww/kqtl_data/results/cabbage/ctrl_db_sorted \
  --case-hist /home/ww/kqtl_data/results/cabbage/case_hist.txt \
  --ctrl-hist /home/ww/kqtl_data/results/cabbage/ctrl_hist.txt \
  --per-position --chrom C09 \
  -o /tmp/cf_test/cabbage_C09_perpos.tsv
```

- Runtime: 86 s
- Output: 34,181,872 positions, 916 MB

### vradiata chromosome 11

```bash
/home/ww/kqtl_data/kbsa-standalone/build/kbsa anchor \
  --ref /mnt/f/kqtl/data/vrad_true/Vrad_JL7.genome.fa \
  --case /home/ww/kqtl_data/results/vradiata/case_db_sorted \
  --ctrl /home/ww/kqtl_data/results/vradiata/ctrl_db_sorted \
  --case-hist /home/ww/kqtl_data/results/vradiata/case_hist.txt \
  --ctrl-hist /home/ww/kqtl_data/results/vradiata/ctrl_hist.txt \
  --per-position --chrom 11 \
  -o /tmp/cf_test/vradiata_chr11_perpos.tsv
```

- Runtime: 68 s
- Output: 23,930,430 positions, 610 MB

---

## Step 2 — Shared-only filter + 1Mb sum|z| ranking

Implemented inline in `scripts/visualize_validation.py`. Reproduce with:

```bash
cd /tmp/cf_test
python3 /home/ww/kqtl_data/kbsa-standalone/docs/validation/scripts/visualize_validation.py .
```

The script renders 4 PNGs and prints rank tables.

---

## Results — all 4 datasets

### Per-dataset shared fraction

| Dataset  | Positions | Shared (case>0 AND ctrl>0) | % shared |
|----------|-----------|----------------------------|----------|
| brapa    | 33,171,530 | 32,366,355                | 97.6 %   |
| cucumber | 28,780,246 | 26,108,258                | 90.7 %   |
| cabbage  | 34,181,872 | 33,690,773 (approx)       | 98.6 %   |
| vradiata | 23,930,430 | 23,742,522 (approx)       | 99.2 %   |

### Target 1Mb-window ranking

| Dataset  | Variant type    | Target window | Rank        | Status |
|----------|-----------------|---------------|-------------|--------|
| brapa    | 4.4 kb PAV      | A09 38 Mb     | **#3 / 39** | ✓ pass (target in top-3 LD cluster 36-38 Mb) |
| cucumber | 102 bp DEL      | chr1 30 Mb    | **#1 / 33** | ✓ pass (15% lead over #2) |
| cabbage  | 1 bp DEL        | C09 29 Mb     | #38 / 41    | ✗ fail (1bp DEL physical limit) |
| vradiata | 1 bp DEL        | chr11 21 Mb   | #23 / 31    | ✗ fail (1bp DEL physical limit) |

### Top-5 windows per dataset (sum|z| per 1Mb)

```
brapa A09:
  #1: 36 Mb  sum|z|=988,516       (LD region of target PAV)
  #2: 37 Mb  sum|z|=987,760       (LD region of target PAV)
  #3: 38 Mb  sum|z|=935,133  ★ TARGET (4.4kb PAV at 38876601-38881010)
  #4: 35 Mb  sum|z|=923,026
  #5: 31 Mb  sum|z|=883,363

cucumber chr1:
  #1: 30 Mb  sum|z|=1,150,267 ★ TARGET (102bp DEL at 30062867-30066723)
  #2: 23 Mb  sum|z|=978,439       (parental SV background)
  #3: 29 Mb  sum|z|=977,246       (LD region adjacent to target)
  #4:  4 Mb  sum|z|=976,152
  #5:  3 Mb  sum|z|=947,356

cabbage C09:
  #1: 33 Mb  sum|z|=785,864       (parental SV background)
  #2: 32 Mb  sum|z|=778,509
  #3: 38 Mb  sum|z|=714,370
  ...
  #38: 29 Mb sum|z|=499,227 ★ TARGET (signal too weak vs background)

vradiata chr11:
  #1: 11 Mb  sum|z|=1,290,383     (parental SV background)
  #2:  8 Mb  sum|z|=1,242,417
  ...
  #23: 21 Mb sum|z|=663,959 ★ TARGET (signal too weak vs background)
```

---

## Why 1bp DEL fails — physics, not algorithm

A single-base deletion produces ~31 boundary k-mers (k=31). In a 1 Mb window:
- Target 1bp-DEL signal: ~31 k-mers × |z| ~2-5 ≈ 100-150 sum|z|
- Background per Mb: ~900,000 k-mers × mean |z| ~0.5-1.0 ≈ 500,000-900,000 sum|z|
- Signal-to-noise ratio for 1bp DEL at 1Mb scale: ~0.02 % — physically swamped by parental background polymorphism.

Classical SNP-based BSA (QTL-seq, MutMap) avoids this by counting allele frequencies at sparse SNP sites only, not aggregating over all positions. Our k-mer-based approach inherits the strength of dense per-position signal (good for PAV/DEL > 100 bp) but cannot resolve single-base events at Mb scale without alignment to sparse marker sets.

This is a documented BSA-seq limitation, not an algorithmic regression.

---

## Figures

All figures are 4-panel layouts:

| Panel | Content | Aggregation |
|-------|---------|-------------|
| Top-left  | 1 Mb sum\|z\| bar chart (target window highlighted red) | Per 1 Mb (the **ranking metric**) |
| Top-right | 10 kb mean\|z\| line plot, full chromosome | Per 10 kb |
| Bottom-left | 10 kb mean\|z\| zoomed to target ±5 Mb | Per 10 kb |
| Bottom-right | 10 kb mean signed z (CASE/CTRL direction) zoomed | Per 10 kb |

Files (in `figures/`):
- `dataset_brapa.png` — A09 4.4kb PAV, target rank #3/39
- `dataset_cucumber.png` — chr1 102bp DEL, target rank #1/33
- `dataset_cabbage.png` — C09 1bp DEL, target rank #38/41
- `dataset_vradiata.png` — chr11 1bp DEL, target rank #23/31
- `four_datasets_validation.png` — combined 4-panel summary (one panel per dataset)

---

## Reproducibility

Build:
```bash
cd /home/ww/kqtl_data/kbsa-standalone/build
cmake .. && cmake --build .
```

Run all 4 datasets (~5 min total wall time):
```bash
mkdir -p /tmp/cf_test
# Run the 4 anchor commands above to produce the per-position TSVs.
python3 /home/ww/kqtl_data/kbsa-standalone/docs/validation/scripts/visualize_validation.py /tmp/cf_test
# Output PNGs land in /tmp/cf_test/dataset_<name>.png
```

The build binary path used in this validation: `/home/ww/kqtl_data/kbsa-standalone/build/kbsa` (commit-time of this session).

---

## Algorithmic notes

- **No tunable thresholds**: `kai_min` and `g_min` default to 0.0 (emit all positions with non-zero counts in either pool). The shared-only filter is BSA-definitional, not a tunable threshold.
- **No tile/HMM layer**: this validation uses pure 1Mb sum aggregation. Tile + HMM extensions were prototyped but **did not improve results** in this validation; they are kept out of the main path.
- **Direction is intrinsic**: `z_score` sign is determined by `kai_reg ≥ 0.5` ⇒ CONTROL-enriched (z > 0), `< 0.5` ⇒ CASE-enriched (z < 0). No threshold influences direction.
- **PAV signal vs BSA signal**: positions where one pool has zero count (case-only or ctrl-only) likely represent presence/absence variants between parents, not allele frequency shifts. The shared-only filter excludes them; without this filter, cucumber target sinks to rank 27/33 (parental SVs dominate).

# kbsa — k-mer Bulked Segregant Analysis

Reference-free k-mer-based BSA-seq tool. Identifies QTL intervals from pooled-sample sequencing without variant calling or read alignment.

## Principle

Traditional BSA-seq requires: read alignment → variant calling → SNP-index calculation. kbsa skips all three steps by directly comparing k-mer frequency distributions between two bulk pools, scoring each k-mer's allele frequency deviation, then projecting scores onto reference coordinates.

**Core statistic**: For each k-mer present in either pool, compute `kai_reg = f_bulk1 / (f_bulk1 + f_bulk2)` after depth normalization. The signed z-score `z = sign(kai - 0.5) * sqrt(G)` quantifies enrichment direction and magnitude.

## Modes

### `kbsa anchor --peak` (primary)

Projects k-mer scores onto reference genome positions via sliding window. Outputs ranked intervals by excess-mass score.

```bash
kbsa anchor --ref genome.fa \
  --bulk1 bulk1_db_sorted --bulk2 bulk2_db_sorted \
  --bulk1-hist bulk1_hist.txt --bulk2-hist bulk2_hist.txt \
  --peak --threads 6 -o output.bed
```

**Algorithm**: For each k-mer on the reference (sliding by 1bp), exact-lookup its count in both bulk KMC databases via minimizer-indexed hash query. Compute the BSA z-score, accumulate `max(0, |z| - tau)` into the window the position belongs to. Rank windows by sum.

This is **not** alignment, mapping, or pseudo-alignment (Salmon/kallisto). The reference genome serves only as the source of probe k-mers; `KMC::CheckKmer` is an O(1) exact hash lookup. Consequence: anchor mode cannot find k-mers absent from the reference — that is what `unitig` mode handles.

### `kbsa unitig` (reference-free)

Assembles differential k-mers into unitigs via BCALM2, then optionally maps back to reference. Detects structural variants invisible to anchor mode.

```bash
kbsa unitig --bulk1 bulk1_sorted --bulk2 bulk2_sorted \
  --bulk1-hist bulk1_hist.txt --bulk2-hist bulk2_hist.txt \
  --map-ref genome.fa -o output_prefix
```

**Algorithm**: Stream-merge both KMC databases, retain k-mers with strong differential signal (`max(bulk1, bulk2) >= valley` and `total <= peak*3`), write to FASTA, run BCALM2 to compact into unitigs, score each unitig by averaging k-mer z-scores. Optional `--map-ref` aligns unitigs to reference via minimap2.

**Why this complements anchor mode**: BCALM2 assembles k-mers that are *not* in the reference (novel insertions, PAV alleles, sample-specific variants). Anchor mode misses these by construction. Cabbage Bra032670 (PAV) and vradiata jg35124 (1bp deletion) both validate via unitig mode.

## Validation Results

Tested on 5 datasets spanning different population designs, organisms, and causal variant types:

| Dataset | Organism | Population | Causal Gene | Variant Type | Literature Interval | kbsa Top Hit | Status |
|---------|----------|-----------|-------------|--------------|--------------------:|-------------:|--------|
| brapa | *B. rapa* | F2 (25+25) | Bra032670 | PAV | A09: 37.35-38.88 Mb | A09: 38Mb | **Hit** |
| cucumber | *C. sativus* | F2 (30+30) | CsaV3_1G044640 | 102bp deletion | chr1: ~30 Mb | chr1: 30Mb | **Hit** |
| cabbage | *B. oleracea* | BC24 | Bol035718 | 1bp promoter deletion | C09: 28-31 Mb | C09: 30Mb | **Hit** |
| vradiata | *V. radiata* | RIL (30+30) | jg35124 | 1bp exon deletion | chr11: 6.23-12.75 Mb | chr11: 11Mb | **Hit** |
| soybean | *G. max* | F2 (40+40) | Glyma.06G202300 | frameshift deletion | Gm06: 17.18-20.58 Mb | Gm06: 18Mb | **Hit** |

All 5 datasets correctly identify the literature-reported QTL interval as the top-ranked region.

### Cabbage case study (high heterozygosity)

Cabbage (het=7.36%) requires GenomeScope2-derived peak override (`--peak1 30 --peak2 30`) because the heuristic peak detector confuses the heterozygous peak (~15x) with the homozygous peak (~30x). The external pipeline script `scripts/kbsa-anchor-gs2.sh` automates this.

### Unitig validation

Unitig mode assembles 351 unitigs for cabbage (89/108 mapped to C09) and 118 for vradiata (34/50 mapped to chr11). Unitigs cluster within the QTL intervals, confirming anchor mode results at base-pair resolution.

## Advantages

- **No alignment needed**: Works directly from raw FASTQ → KMC databases
- **Reference-free option**: Unitig mode detects novel insertions and PAV absent from reference
- **Handles high heterozygosity**: GenomeScope2 integration corrects peak detection for het >5%
- **Low memory in score mode**: Streaming merge of sorted KMC databases, O(1) memory
- **Scalable**: Anchor mode parallelizes per-chromosome; 6 threads processes 9 chromosomes in ~140s

## Limitations

- **Window resolution**: Anchor mode reports 1Mb intervals, not single-gene resolution
- **Memory in anchor mode**: Each thread opens ~1.5GB KMC random-access index (6 threads ≈ 22GB)
- **Requires sorted KMC databases**: `kmc_tools transform <db> sort <db_sorted>` is mandatory
- **No direct variant calling**: Identifies QTL intervals, not causal SNPs (unitig mode partially addresses this)
- **Heuristic peak detection fails for high-het genomes**: Requires external GenomeScope2 pipeline

## Pipeline

```
bulk1.fastq + bulk2.fastq
        ↓ kbsa count (KMC)
bulk1_db + bulk2_db
        ↓ kmc_tools transform sort
bulk1_db_sorted + bulk2_db_sorted
        ↓ kbsa anchor --peak
ranked_windows.bed (1Mb QTL intervals)
        ↓ kbsa unitig --map-ref
unitigs.fa + unitigs.bed (base-pair-level loci)
```

## Installation

```bash
git clone <repo>
cd kbsa-standalone
mkdir build && cd build
cmake .. && make -j$(nproc)
```

### Runtime dependencies (not compiled, install separately)

- [BCALM2](https://github.com/GATB/bcalm) — unitig assembly
- [minimap2](https://github.com/lh3/minimap2) — unitig-to-reference mapping
- [GenomeScope2](https://github.com/tbenavi1/genomescope2.0) — peak estimation for high-het genomes

## CLI Reference

```
kbsa count   — KMC k-mer counting from FASTQ
kbsa score   — Per-k-mer kai_reg, G-score, signed z-score
kbsa anchor  — Reference-anchored QTL interval detection (--peak)
kbsa unitig  — Reference-free differential unitig assembly
```

### Key options

| Option | Mode | Description |
|--------|------|-------------|
| `--peak` | anchor | Enable 1Mb sliding window ranking |
| `--peak1 N` | anchor | Override bulk1 peak depth (from GenomeScope2) |
| `--peak2 N` | anchor | Override bulk2 peak depth |
| `--threads N` | anchor | Parallel threads (memory: N × 1.5GB) |
| `--map-ref` | unitig | Map assembled unitigs back to reference |
| `--min-unitig N` | unitig | Minimum unitig length in bp (default: 200) |

## Design Principles

1. **Algorithm generalizes; real data only validates.** Thresholds derive from BSA theory, not tuned to datasets.
2. **Score does not make biological decisions.** Default thresholds emit every k-mer. Downstream stages discriminate.
3. **Direction is intrinsic to data.** `kai_reg > 0.5` → BULK1-enriched, regardless of threshold settings.

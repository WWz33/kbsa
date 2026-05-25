#!/usr/bin/env python3
"""
Simulate F2 BSA-seq experiment from scratch.

Generates:
  - Synthetic reference chromosome with realistic GC content
  - Parental SNP variants (random background + 1 causal)
  - F2 population with recombination
  - Two bulk pools selected by causal locus genotype
  - Per-haplotype FASTA files ready for read simulation
"""

import argparse
import os
import random
import sys
from pathlib import Path


def generate_reference(length: int, gc_content: float, seed: int) -> str:
    """Generate synthetic chromosome with target GC content."""
    rng = random.Random(seed)
    gc_prob = gc_content / 2
    at_prob = (1 - gc_content) / 2
    weights = [at_prob, gc_prob, gc_prob, at_prob]  # A, C, G, T
    bases = "ACGT"
    return "".join(rng.choices(bases, weights=weights, k=length))


def introduce_snps(ref_seq: str, n_snps: int, causal_pos: int, seed: int):
    """
    Place random background SNPs + 1 causal SNP.
    Returns: list of (pos, ref_base, alt_base) sorted by position.
    """
    rng = random.Random(seed)
    length = len(ref_seq)
    bases = "ACGT"

    positions = set()
    positions.add(causal_pos)

    while len(positions) < n_snps + 1:
        p = rng.randint(0, length - 1)
        positions.add(p)

    variants = []
    for pos in sorted(positions):
        ref_base = ref_seq[pos]
        alt_choices = [b for b in bases if b != ref_base]
        alt_base = rng.choice(alt_choices)
        is_causal = (pos == causal_pos)
        variants.append((pos, ref_base, alt_base, is_causal))

    return variants


def simulate_f2_gamete(variant_positions, chrom_length, recomb_rate_per_bp, rng):
    """
    Simulate one F1 gamete (meiosis with crossovers).
    F1 is heterozygous at all variant positions (0|1).
    Returns: list of alleles (0 or 1) for each variant position.
    """
    # Draw crossover positions (Poisson process)
    expected_crossovers = chrom_length * recomb_rate_per_bp
    n_crossovers = rng.poisson(expected_crossovers)
    crossover_positions = sorted(rng.integers(0, chrom_length, size=n_crossovers))

    # Start on random haplotype
    current_hap = rng.integers(0, 2)
    alleles = []
    co_idx = 0

    for pos in variant_positions:
        while co_idx < len(crossover_positions) and crossover_positions[co_idx] <= pos:
            current_hap = 1 - current_hap
            co_idx += 1
        alleles.append(current_hap)

    return alleles


def simulate_f2_population(variants, chrom_length, recomb_rate, n_f2, seed):
    """
    Simulate n_f2 F2 individuals.
    Each F2 = two independent F1 gametes.
    Returns: list of (hap1_alleles, hap2_alleles) per individual.
    """
    import numpy as np
    rng = np.random.default_rng(seed)

    variant_positions = [v[0] for v in variants]
    population = []

    for _ in range(n_f2):
        hap1 = simulate_f2_gamete(variant_positions, chrom_length, recomb_rate, rng)
        hap2 = simulate_f2_gamete(variant_positions, chrom_length, recomb_rate, rng)
        population.append((hap1, hap2))

    return population


def select_bulks(population, variants, causal_idx, bulk_size):
    """
    Select individuals for two bulks based on causal locus genotype.
    Bulk1: homozygous REF (0/0) at causal locus
    Bulk2: homozygous ALT (1/1) at causal locus
    """
    bulk1_candidates = []
    bulk2_candidates = []

    for i, (hap1, hap2) in enumerate(population):
        geno = hap1[causal_idx] + hap2[causal_idx]
        if geno == 0:
            bulk1_candidates.append(i)
        elif geno == 2:
            bulk2_candidates.append(i)

    if len(bulk1_candidates) < bulk_size:
        print(f"WARNING: only {len(bulk1_candidates)} homREF individuals "
              f"(need {bulk_size}). Using all.", file=sys.stderr)
        bulk1 = bulk1_candidates
    else:
        bulk1 = bulk1_candidates[:bulk_size]

    if len(bulk2_candidates) < bulk_size:
        print(f"WARNING: only {len(bulk2_candidates)} homALT individuals "
              f"(need {bulk_size}). Using all.", file=sys.stderr)
        bulk2 = bulk2_candidates
    else:
        bulk2 = bulk2_candidates[:bulk_size]

    return bulk1, bulk2


def apply_variants_to_seq(ref_seq: str, variants, alleles) -> str:
    """Apply variant alleles to reference sequence."""
    seq = list(ref_seq)
    for i, (pos, ref_base, alt_base, _) in enumerate(variants):
        if alleles[i] == 1:
            seq[pos] = alt_base
    return "".join(seq)


def write_fasta(filepath: str, seq_name: str, sequence: str, line_width=80):
    """Write sequence to FASTA file."""
    with open(filepath, "w") as f:
        f.write(f">{seq_name}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Simulate F2 BSA bulks")
    parser.add_argument("--chrom-length", type=int, default=26_000_000,
                        help="Chromosome length in bp (default: 26Mb)")
    parser.add_argument("--gc-content", type=float, default=0.37,
                        help="GC content (default: 0.37)")
    parser.add_argument("--snp-density", type=float, default=1/500,
                        help="Background SNP density per bp (default: 1/500)")
    parser.add_argument("--causal-pos", type=int, default=None,
                        help="Causal SNP position (default: middle of chrom)")
    parser.add_argument("--recomb-rate", type=float, default=1.5e-8,
                        help="Recombination rate per bp (default: 1.5e-8, ~1.5cM/Mb)")
    parser.add_argument("--n-f2", type=int, default=2000,
                        help="Total F2 population size (default: 2000)")
    parser.add_argument("--bulk-size", type=int, default=40,
                        help="Individuals per bulk (default: 40)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Output directory")
    args = parser.parse_args()

    if args.causal_pos is None:
        args.causal_pos = args.chrom_length // 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[sim] Generating {args.chrom_length/1e6:.1f}Mb reference (GC={args.gc_content})",
          file=sys.stderr)
    ref_seq = generate_reference(args.chrom_length, args.gc_content, args.seed)

    n_background_snps = int(args.chrom_length * args.snp_density)
    print(f"[sim] Introducing {n_background_snps} background SNPs + 1 causal @ {args.causal_pos}",
          file=sys.stderr)
    variants = introduce_snps(ref_seq, n_background_snps, args.causal_pos, args.seed + 1)

    causal_idx = next(i for i, v in enumerate(variants) if v[3])
    print(f"[sim] Causal SNP: pos={variants[causal_idx][0]} "
          f"{variants[causal_idx][1]}>{variants[causal_idx][2]}", file=sys.stderr)

    print(f"[sim] Simulating {args.n_f2} F2 individuals "
          f"(recomb={args.recomb_rate*1e8:.1f} cM/Mb)", file=sys.stderr)
    population = simulate_f2_population(
        variants, args.chrom_length, args.recomb_rate, args.n_f2, args.seed + 2)

    print(f"[sim] Selecting bulks (size={args.bulk_size})", file=sys.stderr)
    bulk1_idx, bulk2_idx = select_bulks(
        population, variants, causal_idx, args.bulk_size)
    print(f"[sim] Bulk1 (homREF): {len(bulk1_idx)} individuals", file=sys.stderr)
    print(f"[sim] Bulk2 (homALT): {len(bulk2_idx)} individuals", file=sys.stderr)

    # Write reference
    ref_path = outdir / "ref.fa"
    write_fasta(str(ref_path), "chr1", ref_seq)
    print(f"[sim] Reference written: {ref_path}", file=sys.stderr)

    # Write variant truth table
    truth_path = outdir / "variants_truth.tsv"
    with open(truth_path, "w") as f:
        f.write("pos\tref\talt\tis_causal\n")
        for pos, ref_b, alt_b, is_causal in variants:
            f.write(f"{pos}\t{ref_b}\t{alt_b}\t{int(is_causal)}\n")

    # Write haplotype FASTAs for each bulk
    for bulk_name, indices in [("bulk1", bulk1_idx), ("bulk2", bulk2_idx)]:
        bulk_dir = outdir / "haps" / bulk_name
        bulk_dir.mkdir(parents=True, exist_ok=True)

        for ind_idx in indices:
            hap1_alleles, hap2_alleles = population[ind_idx]
            # Haplotype 1
            seq1 = apply_variants_to_seq(ref_seq, variants, hap1_alleles)
            write_fasta(str(bulk_dir / f"ind{ind_idx:04d}.h1.fa"),
                        f"chr1_ind{ind_idx}_h1", seq1)
            # Haplotype 2
            seq2 = apply_variants_to_seq(ref_seq, variants, hap2_alleles)
            write_fasta(str(bulk_dir / f"ind{ind_idx:04d}.h2.fa"),
                        f"chr1_ind{ind_idx}_h2", seq2)

    print(f"[sim] Haplotype FASTAs written to {outdir}/haps/", file=sys.stderr)
    print(f"[sim] Done. Causal position: {args.causal_pos}", file=sys.stderr)

    # Write metadata
    meta_path = outdir / "sim_metadata.txt"
    with open(meta_path, "w") as f:
        f.write(f"chrom_length={args.chrom_length}\n")
        f.write(f"gc_content={args.gc_content}\n")
        f.write(f"n_background_snps={n_background_snps}\n")
        f.write(f"causal_pos={args.causal_pos}\n")
        f.write(f"causal_ref={variants[causal_idx][1]}\n")
        f.write(f"causal_alt={variants[causal_idx][2]}\n")
        f.write(f"recomb_rate={args.recomb_rate}\n")
        f.write(f"n_f2={args.n_f2}\n")
        f.write(f"bulk_size={args.bulk_size}\n")
        f.write(f"bulk1_n={len(bulk1_idx)}\n")
        f.write(f"bulk2_n={len(bulk2_idx)}\n")
        f.write(f"seed={args.seed}\n")


if __name__ == "__main__":
    main()

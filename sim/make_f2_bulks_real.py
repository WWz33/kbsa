#!/usr/bin/env python3
"""
Simulate F2 BSA-seq experiment using a REAL reference chromosome.

Key difference from synthetic version: uses actual genomic sequence,
preserving GC structure, repeats, and k-mer complexity.
"""

import argparse
import os
import random
import sys
from pathlib import Path


def load_fasta(filepath: str) -> str:
    """Load single-sequence FASTA into string."""
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    return "".join(seq_parts)


def generate_random_seq(length: int, rng) -> str:
    """Generate random DNA sequence."""
    return "".join(rng.choice("ACGT") for _ in range(length))


def make_causal_variant(ref_seq: str, causal_pos: int, causal_type: str,
                        causal_size: int, rng):
    """
    Create causal variant of specified type.
    Returns: (pos, ref_allele, alt_allele)
    """
    if causal_type == "snp":
        ref_base = ref_seq[causal_pos]
        alt_base = rng.choice([b for b in "ACGT" if b != ref_base])
        return (causal_pos, ref_base, alt_base)

    elif causal_type == "insertion":
        # Insert novel sequence after causal_pos (ref_allele = "", alt = novel seq)
        ref_allele = ref_seq[causal_pos]
        alt_allele = ref_allele + generate_random_seq(causal_size, rng)
        return (causal_pos, ref_allele, alt_allele)

    elif causal_type == "deletion":
        # Delete causal_size bp starting at causal_pos
        ref_allele = ref_seq[causal_pos:causal_pos + causal_size]
        alt_allele = ref_seq[causal_pos]  # keep anchor base
        return (causal_pos, ref_allele, alt_allele)

    elif causal_type == "pav":
        # Presence/absence: large deletion in ALT (gene-scale)
        ref_allele = ref_seq[causal_pos:causal_pos + causal_size]
        alt_allele = ref_seq[causal_pos]
        return (causal_pos, ref_allele, alt_allele)

    else:
        raise ValueError(f"Unknown causal_type: {causal_type}")


def introduce_variants(ref_seq: str, n_snps: int, causal_pos: int,
                       causal_type: str, causal_size: int, seed: int):
    """
    Place random background SNPs + 1 causal variant (SNP/indel/insertion/PAV).
    Returns: list of (pos, ref_allele, alt_allele, is_causal) sorted by position.
    """
    rng = random.Random(seed)
    length = len(ref_seq)
    bases = "ACGT"

    valid_positions = []
    for p in range(length):
        if ref_seq[p] in bases:
            valid_positions.append(p)

    # Ensure causal position is valid
    if ref_seq[causal_pos] not in bases:
        for offset in range(1, 1000):
            if causal_pos + offset < length and ref_seq[causal_pos + offset] in bases:
                causal_pos = causal_pos + offset
                break
            if causal_pos - offset >= 0 and ref_seq[causal_pos - offset] in bases:
                causal_pos = causal_pos - offset
                break

    # Determine exclusion zone around causal variant (avoid overlapping variants)
    if causal_type == "snp":
        causal_span = 1
    else:
        causal_span = max(causal_size, 1)
    exclusion = set(range(causal_pos, causal_pos + causal_span))

    # Sample background SNP positions (excluding causal zone)
    valid_set = set(valid_positions) - exclusion
    valid_list = list(valid_set)

    if n_snps > len(valid_list):
        n_snps = len(valid_list)
        print(f"WARNING: reduced to {n_snps} background SNPs", file=sys.stderr)

    sampled = sorted(rng.sample(valid_list, n_snps))

    # Build variant list: background SNPs
    variants = []
    for pos in sampled:
        ref_base = ref_seq[pos]
        alt_choices = [b for b in bases if b != ref_base]
        alt_base = rng.choice(alt_choices)
        variants.append((pos, ref_base, alt_base, False))

    # Add causal variant
    causal_ref, causal_alt = make_causal_variant(
        ref_seq, causal_pos, causal_type, causal_size, rng)[1:]
    causal_entry = (causal_pos, causal_ref, causal_alt, True)

    # Insert causal in sorted order
    import bisect
    positions = [v[0] for v in variants]
    idx = bisect.bisect_left(positions, causal_pos)
    variants.insert(idx, causal_entry)

    return variants, causal_pos


def simulate_f2_gamete(variant_positions, chrom_length, recomb_rate_per_bp, rng):
    """
    Simulate one F1 gamete (meiosis with crossovers).
    Returns: list of alleles (0 or 1) for each variant position.
    """
    expected_crossovers = chrom_length * recomb_rate_per_bp
    n_crossovers = rng.poisson(expected_crossovers)
    crossover_positions = sorted(rng.integers(0, chrom_length, size=n_crossovers))

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
    """Simulate n_f2 F2 individuals."""
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
    """Select individuals for two bulks based on causal locus genotype."""
    bulk1_candidates = []
    bulk2_candidates = []

    for i, (hap1, hap2) in enumerate(population):
        geno = hap1[causal_idx] + hap2[causal_idx]
        if geno == 0:
            bulk1_candidates.append(i)
        elif geno == 2:
            bulk2_candidates.append(i)

    if len(bulk1_candidates) < bulk_size:
        print(f"WARNING: only {len(bulk1_candidates)} homREF (need {bulk_size})",
              file=sys.stderr)
        bulk1 = bulk1_candidates
    else:
        bulk1 = bulk1_candidates[:bulk_size]

    if len(bulk2_candidates) < bulk_size:
        print(f"WARNING: only {len(bulk2_candidates)} homALT (need {bulk_size})",
              file=sys.stderr)
        bulk2 = bulk2_candidates
    else:
        bulk2 = bulk2_candidates[:bulk_size]

    return bulk1, bulk2


def apply_variants_to_seq(ref_seq: str, variants, alleles) -> str:
    """Apply variant alleles to reference sequence (supports indels). O(n) forward scan."""
    parts = []
    prev_end = 0
    for i in range(len(variants)):
        if alleles[i] == 1:
            pos, ref_allele, alt_allele, _ = variants[i]
            parts.append(ref_seq[prev_end:pos])
            parts.append(alt_allele)
            prev_end = pos + len(ref_allele)
    parts.append(ref_seq[prev_end:])
    return "".join(parts)


def write_fasta(filepath: str, seq_name: str, sequence: str, line_width=80):
    """Write sequence to FASTA file."""
    with open(filepath, "w") as f:
        f.write(f">{seq_name}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Simulate F2 BSA on real reference")
    parser.add_argument("--ref", type=str, required=True,
                        help="Reference FASTA (single chromosome)")
    parser.add_argument("--chrom-name", type=str, default="chr1",
                        help="Chromosome name for output")
    parser.add_argument("--snp-density", type=float, default=1/500,
                        help="Background SNP density per bp (default: 1/500)")
    parser.add_argument("--causal-pos", type=int, required=True,
                        help="Causal variant position (0-based)")
    parser.add_argument("--causal-type", type=str, default="snp",
                        choices=["snp", "insertion", "deletion", "pav"],
                        help="Causal variant type (default: snp)")
    parser.add_argument("--causal-size", type=int, default=1,
                        help="Size of causal variant in bp (ignored for snp)")
    parser.add_argument("--recomb-rate", type=float, default=1.5e-8,
                        help="Recombination rate per bp (default: 1.5e-8)")
    parser.add_argument("--n-f2", type=int, default=2000,
                        help="Total F2 population size")
    parser.add_argument("--bulk-size", type=int, default=40,
                        help="Individuals per bulk")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[sim] Loading reference: {args.ref}", file=sys.stderr)
    ref_seq = load_fasta(args.ref)
    chrom_length = len(ref_seq)
    print(f"[sim] Chromosome length: {chrom_length/1e6:.1f}Mb", file=sys.stderr)

    n_background_snps = int(chrom_length * args.snp_density)
    print(f"[sim] Introducing {n_background_snps} background SNPs + 1 causal "
          f"{args.causal_type}({args.causal_size}bp) @ {args.causal_pos}",
          file=sys.stderr)
    variants, causal_pos = introduce_variants(
        ref_seq, n_background_snps, args.causal_pos,
        args.causal_type, args.causal_size, args.seed + 1)

    causal_idx = next(i for i, v in enumerate(variants) if v[3])
    causal_v = variants[causal_idx]
    ref_display = causal_v[1] if len(causal_v[1]) <= 10 else f"{causal_v[1][:5]}..({len(causal_v[1])}bp)"
    alt_display = causal_v[2] if len(causal_v[2]) <= 10 else f"{causal_v[2][:5]}..({len(causal_v[2])}bp)"
    print(f"[sim] Causal {args.causal_type}: pos={causal_v[0]} "
          f"{ref_display}>{alt_display}", file=sys.stderr)

    print(f"[sim] Simulating {args.n_f2} F2 individuals "
          f"(recomb={args.recomb_rate*1e8:.1f} cM/Mb)", file=sys.stderr)
    population = simulate_f2_population(
        variants, chrom_length, args.recomb_rate, args.n_f2, args.seed + 2)

    print(f"[sim] Selecting bulks (size={args.bulk_size})", file=sys.stderr)
    bulk1_idx, bulk2_idx = select_bulks(
        population, variants, causal_idx, args.bulk_size)
    print(f"[sim] Bulk1 (homREF): {len(bulk1_idx)} individuals", file=sys.stderr)
    print(f"[sim] Bulk2 (homALT): {len(bulk2_idx)} individuals", file=sys.stderr)

    # Write reference (just this chromosome)
    ref_path = outdir / "ref.fa"
    write_fasta(str(ref_path), args.chrom_name, ref_seq)
    print(f"[sim] Reference written: {ref_path}", file=sys.stderr)

    # Write variant truth table
    truth_path = outdir / "variants_truth.tsv"
    with open(truth_path, "w") as f:
        f.write("pos\tref\talt\tis_causal\n")
        for pos, ref_b, alt_b, is_causal in variants:
            f.write(f"{pos}\t{ref_b}\t{alt_b}\t{int(is_causal)}\n")

    # Write haplotype FASTAs
    for bulk_name, indices in [("bulk1", bulk1_idx), ("bulk2", bulk2_idx)]:
        bulk_dir = outdir / "haps" / bulk_name
        bulk_dir.mkdir(parents=True, exist_ok=True)

        for ind_idx in indices:
            hap1_alleles, hap2_alleles = population[ind_idx]
            seq1 = apply_variants_to_seq(ref_seq, variants, hap1_alleles)
            write_fasta(str(bulk_dir / f"ind{ind_idx:04d}.h1.fa"),
                        f"{args.chrom_name}_ind{ind_idx}_h1", seq1)
            seq2 = apply_variants_to_seq(ref_seq, variants, hap2_alleles)
            write_fasta(str(bulk_dir / f"ind{ind_idx:04d}.h2.fa"),
                        f"{args.chrom_name}_ind{ind_idx}_h2", seq2)

    print(f"[sim] Haplotype FASTAs written to {outdir}/haps/", file=sys.stderr)

    # Write metadata
    meta_path = outdir / "sim_metadata.txt"
    with open(meta_path, "w") as f:
        f.write(f"ref_source={args.ref}\n")
        f.write(f"chrom_name={args.chrom_name}\n")
        f.write(f"chrom_length={chrom_length}\n")
        f.write(f"n_background_snps={n_background_snps}\n")
        f.write(f"causal_pos={causal_pos}\n")
        f.write(f"causal_type={args.causal_type}\n")
        f.write(f"causal_size={args.causal_size}\n")
        f.write(f"causal_ref={causal_v[1]}\n")
        f.write(f"causal_alt={causal_v[2]}\n")
        f.write(f"recomb_rate={args.recomb_rate}\n")
        f.write(f"n_f2={args.n_f2}\n")
        f.write(f"bulk_size={args.bulk_size}\n")
        f.write(f"bulk1_n={len(bulk1_idx)}\n")
        f.write(f"bulk2_n={len(bulk2_idx)}\n")
        f.write(f"seed={args.seed}\n")

    print(f"[sim] Done. Causal position: {causal_pos}", file=sys.stderr)


if __name__ == "__main__":
    main()

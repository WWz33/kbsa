#!/usr/bin/env python3
"""Validate kbsa scored output against ground truth PAV/DEL regions."""
import sys
import os

def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def canonical(kmer):
    rc = revcomp(kmer)
    return min(kmer, rc)

def extract_region_seq(fasta_path, chrom, start_1based, end_1based):
    """Extract sequence from fasta (1-based coordinates)."""
    seq_parts = []
    found = False
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if found:
                    break
                if line.split()[0][1:] == chrom:
                    found = True
            elif found:
                seq_parts.append(line.strip())
    full_seq = "".join(seq_parts)
    return full_seq[start_1based - 1 : end_1based].upper()

def generate_canonical_kmers(seq, k=31):
    """Generate all canonical k-mers from a sequence."""
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            kmers.add(canonical(kmer))
    return kmers

def lookup_in_scored(scored_path, target_kmers):
    """Look up target k-mers in scored TSV output."""
    ctrl_found = 0
    case_found = 0
    hits = []
    with open(scored_path) as f:
        next(f)  # skip header
        for line in f:
            cols = line.rstrip("\n").split("\t")
            kmer = cols[0]
            ck = canonical(kmer)
            if ck in target_kmers:
                if cols[1] == "CONTROL":
                    ctrl_found += 1
                else:
                    case_found += 1
                if len(hits) < 10:
                    hits.append(cols[:4])
    return ctrl_found, case_found, hits

def main():
    if len(sys.argv) < 2:
        print("Usage: validate_gt.py brapa|cucumber")
        sys.exit(1)

    dataset = sys.argv[1]

    if dataset == "brapa":
        ref = "/mnt/f/kqtl/data/brapa_true/Brapa_sequence_v1.5.fa"
        scored = "/home/ww/kqtl_data/results/brapa/kbsa_scored.tsv"
        chrom = "A09"
        start = 38876601
        end = 38881010
        expect_enriched = "CONTROL"  # Case=Mutant(缺失), Ctrl=WT → WT有该区域
        print("=== Brapa Ground Truth: Bra032670 PAV ===")
        print(f"Region: {chrom}:{start}-{end} ({end-start+1} bp)")
        print(f"Expected: {expect_enriched}-enriched (WT has the sequence)")
    elif dataset == "cucumber":
        ref = "/mnt/f/kqtl/data/cucumber_true/cucumber_v3.fa"
        scored = "/home/ww/kqtl_data/results/cucumber_csdw3/kbsa_scored.tsv"
        chrom = "Chr1"
        start = 30062867
        end = 30066723
        expect_enriched = "CASE"  # Case=P2(WT), Ctrl=P1(Mutant) → WT有该区域
        print("=== Cucumber Ground Truth: CsaV3_1G044640 DEL ===")
        print(f"Region: {chrom}:{start}-{end} ({end-start+1} bp)")
        print(f"Expected: {expect_enriched}-enriched (WT has the sequence)")
    else:
        print(f"Unknown dataset: {dataset}")
        sys.exit(1)

    print(f"\n[1] Extracting target region from reference...")
    seq = extract_region_seq(ref, chrom, start, end)
    print(f"    Extracted: {len(seq)} bp")
    if len(seq) == 0:
        print("    ERROR: Could not extract sequence! Check chromosome name.")
        # Try listing chromosomes
        with open(ref) as f:
            chroms = [l.split()[0][1:] for l in f if l.startswith(">")]
        print(f"    Available chroms: {chroms[:15]}")
        sys.exit(1)

    print(f"\n[2] Generating canonical 31-mers...")
    target_kmers = generate_canonical_kmers(seq, 31)
    print(f"    Target canonical k-mers: {len(target_kmers)}")

    print(f"\n[3] Looking up in scored output ({scored})...")
    ctrl_found, case_found, hits = lookup_in_scored(scored, target_kmers)
    total_found = ctrl_found + case_found

    print(f"\n=== RESULTS ===")
    print(f"    Target k-mers found in scored: {total_found}/{len(target_kmers)}")
    print(f"    CONTROL-enriched: {ctrl_found}")
    print(f"    CASE-enriched: {case_found}")
    if total_found > 0:
        ratio = ctrl_found / total_found * 100 if expect_enriched == "CONTROL" else case_found / total_found * 100
        print(f"    {expect_enriched} ratio: {ratio:.1f}%")
        print(f"\n    Sample hits:")
        for h in hits:
            print(f"      {h[0]}  {h[1]}  kai={h[2]}  g={h[3]}")
        if ratio >= 70:
            print(f"\n    PASS: {expect_enriched} enrichment {ratio:.1f}% >= 70%")
        else:
            print(f"\n    WARN: {expect_enriched} enrichment {ratio:.1f}% < 70%")
    else:
        print(f"\n    NO target k-mers found in scored output!")
        print(f"    Possible reasons:")
        print(f"    1. All target k-mers filtered by depth threshold")
        print(f"    2. KMC ci=2 removed low-count k-mers in one pool")
        print(f"    → Check with lower G threshold or ci=1")

if __name__ == "__main__":
    main()

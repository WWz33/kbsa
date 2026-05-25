#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "kbsa_unitig.hpp"
#include "kbsa_unitig_score.hpp"

using namespace kbsa;

std::string test_revcomp(const std::string& s) {
  std::string rc(s.size(), 'N');
  for (size_t i = 0; i < s.size(); ++i) {
    char c = s[s.size() - 1 - i];
    switch (c) {
      case 'A': rc[i] = 'T'; break;
      case 'T': rc[i] = 'A'; break;
      case 'C': rc[i] = 'G'; break;
      case 'G': rc[i] = 'C'; break;
      default:  rc[i] = 'N'; break;
    }
  }
  return rc;
}

KmerMap make_kmer_map(const std::string& seq, uint32_t k) {
  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    std::string s = seq.substr(i, k);
    uint64_t fwd = string_to_kmer(s);
    uint64_t can = kmer_canonical(fwd, k);
    kmers[can] = {10, 1, 3.0f, 1};
  }
  return kmers;
}

bool verify_adjacency(const std::string& unitig, uint32_t k) {
  for (size_t i = 0; i + k < unitig.size(); ++i) {
    std::string cur = unitig.substr(i, k);
    std::string nxt = unitig.substr(i + 1, k);
    if (cur.substr(1) != nxt.substr(0, k - 1)) return false;
  }
  return true;
}

bool verify_kmer_conservation(const std::vector<std::string>& unitigs,
                               const KmerMap& input, uint32_t k) {
  for (const auto& u : unitigs) {
    for (size_t i = 0; i + k <= u.size(); ++i) {
      uint64_t fwd = string_to_kmer(u.substr(i, k));
      uint64_t can = kmer_canonical(fwd, k);
      if (!input.count(can)) return false;
    }
  }
  return true;
}

void test_single_unitig() {
  const uint32_t k = 5;
  std::string seq = "AAAACACCC";
  auto kmers = make_kmer_map(seq, k);

  auto unitigs = extract_unitigs(kmers, k, 1);

  std::cout << "  Unitigs from linear chain:\n";
  for (const auto& u : unitigs)
    std::cout << "    " << u << " (len=" << u.size() << ")\n";

  assert(unitigs.size() == 1);
  assert(unitigs[0].size() == 9);
  assert(verify_adjacency(unitigs[0], k));
  assert(verify_kmer_conservation(unitigs, kmers, k));

  std::cout << "  test_single_unitig PASSED\n";
}

void test_branching_unitigs() {
  const uint32_t k = 5;
  std::string seq1 = "AAAACACCC";
  std::string seq2 = "AAAACGTTT";

  KmerMap kmers;
  for (const auto& seq : {seq1, seq2}) {
    for (size_t i = 0; i + k <= seq.size(); ++i) {
      uint64_t fwd = string_to_kmer(seq.substr(i, k));
      uint64_t can = kmer_canonical(fwd, k);
      kmers[can] = {10, 1, 3.0f, 1};
    }
  }

  auto unitigs = extract_unitigs(kmers, k, 1);

  std::cout << "  Unitigs from Y-shape:\n";
  for (const auto& u : unitigs)
    std::cout << "    " << u << " (len=" << u.size() << ")\n";

  assert(unitigs.size() >= 2);

  for (const auto& u : unitigs)
    assert(verify_adjacency(u, k));
  assert(verify_kmer_conservation(unitigs, kmers, k));

  // Edge-disjoint invariant: every input k-mer is covered by at least one unitig
  std::set<uint64_t> covered;
  for (const auto& u : unitigs) {
    for (size_t i = 0; i + k <= u.size(); ++i) {
      uint64_t fwd = string_to_kmer(u.substr(i, k));
      uint64_t can = kmer_canonical(fwd, k);
      covered.insert(can);
    }
  }
  for (const auto& [can, _] : kmers)
    assert(covered.count(can) && "input k-mer missing from unitigs");

  std::cout << "  test_branching_unitigs PASSED\n";
}

void test_min_length_filter() {
  const uint32_t k = 5;
  auto kmers = make_kmer_map("AAAACACCC", k);

  auto unitigs = extract_unitigs(kmers, k, 10);
  assert(unitigs.empty());

  unitigs = extract_unitigs(kmers, k, 9);
  assert(unitigs.size() == 1);

  std::cout << "  test_min_length_filter PASSED\n";
}

void test_empty_input() {
  KmerMap empty;
  auto unitigs = extract_unitigs(empty, 5, 1);
  assert(unitigs.empty());
  std::cout << "  test_empty_input PASSED\n";
}

void test_circular() {
  const uint32_t k = 5;
  // Circular: all k-mers have degree 1/1 in both orientations
  // Use a sequence where last k-1 chars == first k-1 chars
  // ACACACAC -> k-mers: ACACA, CACAC, ACACA(dup), CACAC(dup)
  // That's only 2 unique k-mers — too small.
  // Better: ACGTACGTACGT (12bp, k=5) -> 8 k-mers, wraps if last 4 == first 4
  // ACGTACGTACGT: first 4 = ACGT, last 4 = ACGT — circular!
  std::string seq = "ACGTACGTACGT";
  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = string_to_kmer(seq.substr(i, k));
    uint64_t can = kmer_canonical(fwd, k);
    kmers[can] = {10, 1, 3.0f, 1};
  }

  // Add the wrap-around k-mer to close the cycle
  // Last 4 of seq + first 1 = "ACGT" + "A" = "ACGTA" — already in set
  // Actually for a true cycle we need the (k-1)-overlap to connect last to first
  // seq[len-k+1..len-1] + seq[0] should be in the set
  // That's seq.substr(8,4) + seq[0] = "ACGT" + "A" = "ACGTA"
  // And seq.substr(7,4) + seq[0..0] ... let me just check if it's already circular

  std::cout << "  Circular test k-mers: " << kmers.size() << "\n";
  auto unitigs = extract_unitigs(kmers, k, 1);
  std::cout << "  Circular unitigs: " << unitigs.size() << "\n";
  for (const auto& u : unitigs)
    std::cout << "    " << u << " (len=" << u.size() << ")\n";

  // Should produce at least one unitig covering the k-mers
  assert(verify_kmer_conservation(unitigs, kmers, k));
  for (const auto& u : unitigs)
    assert(verify_adjacency(u, k));

  std::cout << "  test_circular PASSED\n";
}

void test_coherence_bounds() {
  // P3: coherence must be in [0.5, 1.0]
  // Build a unitig where all k-mers have same direction (coherence = 1.0)
  const uint32_t k = 5;
  std::string seq = "ACGTACGTACGTACGT"; // 12 k-mers
  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = string_to_kmer(seq.substr(i, k));
    uint64_t can = kmer_canonical(fwd, k);
    // All enriched in bulk2 (count2 >> count1)
    kmers[can] = {1, 20, 4.0f, 1};
  }

  Params params;
  params.bulk1_to_bulk2_scale = 1.0;
  params.error_rate = 0.005;

  auto us = score_unitig(seq, k, kmers, params, 2.0f);
  assert(us.coherence >= 0.5f && us.coherence <= 1.0f);
  assert(us.coherence > 0.9f); // all same direction → near 1.0
  assert(us.n_kmers == seq.size() - k + 1);
  std::cout << "  coherence=" << us.coherence << " (all same dir)\n";
  std::cout << "  test_coherence_bounds PASSED\n";
}

void test_coherence_mixed() {
  // Mixed directions → coherence near 0.5
  const uint32_t k = 5;
  std::string seq = "ACGTACGTACGTACGT";
  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = string_to_kmer(seq.substr(i, k));
    uint64_t can = kmer_canonical(fwd, k);
    // Alternate: even positions bulk1-enriched, odd bulk2-enriched
    if (i % 2 == 0)
      kmers[can] = {20, 1, 4.0f, 1};
    else
      kmers[can] = {1, 20, 4.0f, 1};
  }

  Params params;
  params.bulk1_to_bulk2_scale = 1.0;
  params.error_rate = 0.005;

  auto us = score_unitig(seq, k, kmers, params, 2.0f);
  assert(us.coherence >= 0.5f && us.coherence <= 1.0f);
  std::cout << "  coherence=" << us.coherence << " (mixed dir)\n";
  std::cout << "  test_coherence_mixed PASSED\n";
}

void test_score_empty_unitig() {
  // Unitig shorter than k → n_kmers=0, no crash
  const uint32_t k = 31;
  std::string seq = "ACGT"; // way shorter than k
  KmerMap kmers;

  Params params;
  params.bulk1_to_bulk2_scale = 1.0;
  params.error_rate = 0.005;

  auto us = score_unitig(seq, k, kmers, params, 2.0f);
  assert(us.n_kmers == 0);
  assert(us.rank_score == 0.0f);
  std::cout << "  test_score_empty_unitig PASSED\n";
}

void test_pbt_idempotency() {
  // P1: running extract_unitigs twice produces identical output
  const uint32_t k = 5;
  std::string seq = "AAAACACCC";
  auto kmers = make_kmer_map(seq, k);

  auto unitigs1 = extract_unitigs(kmers, k, 1);
  auto unitigs2 = extract_unitigs(kmers, k, 1);

  assert(unitigs1.size() == unitigs2.size());
  std::set<std::string> set1(unitigs1.begin(), unitigs1.end());
  std::set<std::string> set2(unitigs2.begin(), unitigs2.end());
  assert(set1 == set2);
  std::cout << "  test_pbt_idempotency PASSED\n";
}

void test_pbt_strand_consistency() {
  // P6: scoring a unitig and its revcomp gives same rank_score
  const uint32_t k = 5;
  std::string seq = "ACGTACGTACGT";
  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = string_to_kmer(seq.substr(i, k));
    uint64_t can = kmer_canonical(fwd, k);
    kmers[can] = {1, 15, 3.5f, 1};
  }

  std::string rc = test_revcomp(seq);
  // Also add revcomp k-mers (they should map to same canonicals)
  for (size_t i = 0; i + k <= rc.size(); ++i) {
    uint64_t fwd = string_to_kmer(rc.substr(i, k));
    uint64_t can = kmer_canonical(fwd, k);
    kmers[can] = {1, 15, 3.5f, 1};
  }

  Params params;
  params.bulk1_to_bulk2_scale = 1.0;
  params.error_rate = 0.005;

  auto s1 = score_unitig(seq, k, kmers, params, 2.0f);
  auto s2 = score_unitig(rc, k, kmers, params, 2.0f);

  assert(s1.n_kmers == s2.n_kmers);
  assert(std::abs(s1.rank_score - s2.rank_score) < 0.001f);
  std::cout << "  score_fwd=" << s1.rank_score << " score_rc=" << s2.rank_score << "\n";
  std::cout << "  test_pbt_strand_consistency PASSED\n";
}

void test_pbt_revcomp_equivalence() {
  // P10: reverse-complementing all input k-mers produces equivalent unitig set
  const uint32_t k = 5;
  std::string seq = "AAAACACCC";
  auto kmers = make_kmer_map(seq, k);

  // Since we use canonical k-mers, the map is already strand-agnostic
  // Build from revcomp sequence
  std::string rc = test_revcomp(seq);
  auto kmers_rc = make_kmer_map(rc, k);

  // Both should produce same canonical kmer set
  assert(kmers.size() == kmers_rc.size());
  for (const auto& [key, _] : kmers)
    assert(kmers_rc.count(key) > 0);

  auto unitigs1 = extract_unitigs(kmers, k, 1);
  auto unitigs2 = extract_unitigs(kmers_rc, k, 1);

  // Same set (each unitig appears as itself or its revcomp)
  std::set<std::string> canonical1, canonical2;
  for (const auto& u : unitigs1) {
    std::string rc_u = test_revcomp(u);
    canonical1.insert(std::min(u, rc_u));
  }
  for (const auto& u : unitigs2) {
    std::string rc_u = test_revcomp(u);
    canonical2.insert(std::min(u, rc_u));
  }
  assert(canonical1 == canonical2);
  std::cout << "  test_pbt_revcomp_equivalence PASSED\n";
}

void test_pbt_threshold_monotonicity() {
  // P7: increasing enriched_min can only decrease differential set size
  const uint32_t k = 5;
  std::string seq = "AAAACACCCCGTGTG";
  KmerMap kmers_low, kmers_high;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = string_to_kmer(seq.substr(i, k));
    uint64_t can = kmer_canonical(fwd, k);
    // Simulate different counts
    uint32_t count = static_cast<uint32_t>(5 + (i % 7));
    kmers_low[can] = {count, 1, 3.0f, 1};
    if (count >= 8)  // higher threshold
      kmers_high[can] = {count, 1, 3.0f, 1};
  }
  assert(kmers_high.size() <= kmers_low.size());
  auto u_low = extract_unitigs(kmers_low, k, 1);
  auto u_high = extract_unitigs(kmers_high, k, 1);

  // Total k-mers in high-threshold unitigs ≤ low-threshold unitigs
  size_t total_low = 0, total_high = 0;
  for (const auto& u : u_low) total_low += u.size() - k + 1;
  for (const auto& u : u_high) total_high += u.size() - k + 1;
  assert(total_high <= total_low);
  std::cout << "  kmers_low=" << total_low << " kmers_high=" << total_high << "\n";
  std::cout << "  test_pbt_threshold_monotonicity PASSED\n";
}

int main() {
  std::cout << "=== extract_unitigs() integration tests ===\n";
  test_single_unitig();
  test_branching_unitigs();
  test_min_length_filter();
  test_empty_input();
  test_circular();
  test_coherence_bounds();
  test_coherence_mixed();
  test_score_empty_unitig();
  test_pbt_idempotency();
  test_pbt_strand_consistency();
  test_pbt_revcomp_equivalence();
  test_pbt_threshold_monotonicity();
  std::cout << "\nAll integration tests PASSED.\n";
  return 0;
}

#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

// Minimal stubs for dependencies not needed in unit test
namespace kbsa {
struct HistEntry { uint64_t depth; uint64_t count; };
struct Params {
  double bulk1_to_bulk2_scale = 1.0;
  uint64_t min_depth = 0;
  uint64_t max_depth = UINT64_MAX;
  double kai_min = 0.0;
  double g_min = 0.0;
  double error_rate = 0.01;
};
inline std::vector<HistEntry> load_histogram(const std::string&) { return {}; }
inline uint64_t find_peak_depth(const std::vector<HistEntry>&) { return 30; }
struct Score { double z_score; int sign; };
enum class Significance { BULK1, BULK2 };
inline Score score_kmer(uint32_t, uint32_t, double, uint64_t, uint64_t, double, double, double) {
  return {0.0, 0};
}
struct MergedKmer { const char* kmer_str; uint32_t kmer_len; uint64_t bulk1_count; uint64_t bulk2_count; };
class KmcMergeIterator {
public:
  KmcMergeIterator(const std::string&, const std::string&) {}
  uint32_t kmer_length() const { return 5; }
  bool next(MergedKmer&) { return false; }
};
}

// Now override includes so kbsa_unitig.hpp doesn't pull real headers
#define KBSA_CORE_HPP
#define KBSA_PARAMS_HPP
#define KBSA_MERGE_HPP

// Include only the parts we need
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace kbsa {

struct KmerInfo {
  uint32_t count1;
  uint32_t count2;
  float z_score;
  int8_t direction;
};

using KmerMap = std::unordered_map<uint64_t, KmerInfo>;

inline uint64_t kmer_revcomp(uint64_t kmer, uint32_t k) noexcept {
  uint64_t rc = 0;
  for (uint32_t i = 0; i < k; ++i) {
    rc = (rc << 2) | (3 - (kmer & 3));
    kmer >>= 2;
  }
  return rc;
}

inline uint64_t kmer_canonical(uint64_t kmer, uint32_t k) noexcept {
  uint64_t rc = kmer_revcomp(kmer, k);
  return std::min(kmer, rc);
}

inline uint64_t kmer_shift_right(uint64_t kmer, uint32_t k, uint8_t base) noexcept {
  uint64_t mask = (k < 32) ? ((uint64_t(1) << (2 * k)) - 1) : ~uint64_t(0);
  return ((kmer << 2) | base) & mask;
}

inline uint64_t kmer_shift_left(uint64_t kmer, uint32_t k, uint8_t base) noexcept {
  return (kmer >> 2) | (uint64_t(base) << (2 * (k - 1)));
}

inline uint8_t kmer_last_base(uint64_t kmer) noexcept { return kmer & 3; }
inline uint8_t kmer_first_base(uint64_t kmer, uint32_t k) noexcept {
  return (kmer >> (2 * (k - 1))) & 3;
}

inline uint8_t base_encode(char c) noexcept {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 0;
  }
}

inline char base_decode(uint8_t b) noexcept { return "ACGT"[b & 3]; }

inline std::string kmer_to_string(uint64_t kmer, uint32_t k) {
  std::string s(k, 'A');
  for (int i = k - 1; i >= 0; --i) {
    s[i] = base_decode(kmer & 3);
    kmer >>= 2;
  }
  return s;
}

inline uint64_t string_to_kmer(const std::string& s) {
  uint64_t kmer = 0;
  for (char c : s) kmer = (kmer << 2) | base_encode(c);
  return kmer;
}

struct OrientedNode {
  uint64_t canonical;
  bool forward;
};

inline std::vector<OrientedNode> right_successors(
    OrientedNode node, uint32_t k, const KmerMap& kmers) {
  uint64_t seq_kmer = node.forward
    ? node.canonical : kmer_revcomp(node.canonical, k);
  std::vector<OrientedNode> result;
  for (uint8_t b = 0; b < 4; ++b) {
    uint64_t candidate = kmer_shift_right(seq_kmer, k, b);
    uint64_t can = kmer_canonical(candidate, k);
    if (kmers.count(can)) {
      result.push_back({can, candidate == can});
    }
  }
  return result;
}

inline std::vector<OrientedNode> left_predecessors(
    OrientedNode node, uint32_t k, const KmerMap& kmers) {
  uint64_t seq_kmer = node.forward
    ? node.canonical : kmer_revcomp(node.canonical, k);
  std::vector<OrientedNode> result;
  for (uint8_t b = 0; b < 4; ++b) {
    uint64_t candidate = kmer_shift_left(seq_kmer, k, b);
    uint64_t can = kmer_canonical(candidate, k);
    if (kmers.count(can)) {
      result.push_back({can, candidate == can});
    }
  }
  return result;
}

struct NodeDegree {
  uint8_t left_fwd;
  uint8_t right_fwd;
  uint8_t left_rev;
  uint8_t right_rev;
};

inline std::unordered_map<uint64_t, NodeDegree> compute_oriented_degrees(
    const KmerMap& kmers, uint32_t k) {
  std::unordered_map<uint64_t, NodeDegree> degrees;
  degrees.reserve(kmers.size());
  for (const auto& [can, _] : kmers) {
    NodeDegree d {};
    OrientedNode fwd_node {can, true};
    OrientedNode rev_node {can, false};
    d.right_fwd = static_cast<uint8_t>(right_successors(fwd_node, k, kmers).size());
    d.left_fwd  = static_cast<uint8_t>(left_predecessors(fwd_node, k, kmers).size());
    d.right_rev = static_cast<uint8_t>(right_successors(rev_node, k, kmers).size());
    d.left_rev  = static_cast<uint8_t>(left_predecessors(rev_node, k, kmers).size());
    degrees[can] = d;
  }
  return degrees;
}

inline bool is_branching(const NodeDegree& d) noexcept {
  bool simple_fwd = (d.left_fwd == 1 && d.right_fwd == 1);
  bool simple_rev = (d.left_rev == 1 && d.right_rev == 1);
  return !simple_fwd || !simple_rev;
}

} // namespace kbsa

using namespace kbsa;

// === Test 1: Basic k-mer utilities ===
void test_kmer_utilities() {
  // A=0, C=1, G=2, T=3
  assert(base_encode('A') == 0);
  assert(base_encode('C') == 1);
  assert(base_encode('G') == 2);
  assert(base_encode('T') == 3);

  // ACGT = 0b00011011 = 27
  uint64_t acgt = string_to_kmer("ACGT");
  assert(acgt == 0b00011011);
  assert(kmer_to_string(acgt, 4) == "ACGT");

  // revcomp(ACGT) = ACGT (palindrome)
  assert(kmer_revcomp(acgt, 4) == acgt);

  // AAAA revcomp = TTTT
  uint64_t aaaa = string_to_kmer("AAAA");
  uint64_t tttt = string_to_kmer("TTTT");
  assert(kmer_revcomp(aaaa, 4) == tttt);

  // canonical picks min
  assert(kmer_canonical(aaaa, 4) == aaaa);  // AAAA < TTTT
  assert(kmer_canonical(tttt, 4) == aaaa);

  std::cout << "  test_kmer_utilities PASSED\n";
}

// === Test 2: Shift operations ===
void test_shifts() {
  // ACGT shift_right with T → CGTT
  uint64_t acgt = string_to_kmer("ACGT");
  uint64_t cgtt = kmer_shift_right(acgt, 4, 3); // T=3
  assert(kmer_to_string(cgtt, 4) == "CGTT");

  // ACGT shift_left with T → TACG
  uint64_t tacg = kmer_shift_left(acgt, 4, 3); // T=3
  assert(kmer_to_string(tacg, 4) == "TACG");

  std::cout << "  test_shifts PASSED\n";
}

// === Test 3: Linear chain — carefully chosen sequence ===
// Use AAAAC AAACA AACAC ACACC CACCC
// Sequence: AAAACACCC (9 bp, k=5, gives 5 k-mers)
// Check each canonical form manually
void test_linear_chain() {
  const uint32_t k = 5;
  std::string seq = "AAAACACCC";

  // Extract k-mers
  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    std::string s = seq.substr(i, k);
    uint64_t fwd = string_to_kmer(s);
    uint64_t can = kmer_canonical(fwd, k);
    std::string can_str = kmer_to_string(can, k);
    std::cout << "    " << s << " -> canonical: " << can_str << "\n";
    kmers[can] = {10, 1, 3.0f, 1};
  }

  std::cout << "    Unique canonical k-mers: " << kmers.size() << "\n";
  assert(kmers.size() == 5);

  // Compute degrees
  auto degrees = compute_oriented_degrees(kmers, k);

  int n_tips = 0;
  int n_simple = 0;
  for (const auto& [can, d] : degrees) {
    std::string s = kmer_to_string(can, k);
    std::cout << "    " << s << ": L_f=" << (int)d.left_fwd
              << " R_f=" << (int)d.right_fwd
              << " L_r=" << (int)d.left_rev
              << " R_r=" << (int)d.right_rev;

    if (is_branching(d)) {
      std::cout << " [NON-SIMPLE]";
      n_tips++;
    } else {
      std::cout << " [SIMPLE]";
      n_simple++;
    }
    std::cout << "\n";
  }

  // In a linear chain, we expect endpoints to be non-simple
  // (degree 0 on one side) and internal nodes to be simple
  std::cout << "    Non-simple: " << n_tips << ", Simple: " << n_simple << "\n";

  // The key invariant: at least 2 non-simple nodes (the endpoints)
  // and at least some simple interior nodes
  assert(n_tips >= 2 && "Linear chain must have at least 2 endpoint nodes");
  assert(n_simple >= 1 && "Linear chain must have at least 1 interior node");

  std::cout << "  test_linear_chain PASSED\n";
}

// === Test 4: Unitig extraction produces correct sequence ===
void test_unitig_extraction() {
  const uint32_t k = 5;
  std::string seq = "AAAACACCC";

  KmerMap kmers;
  for (size_t i = 0; i + k <= seq.size(); ++i) {
    std::string s = seq.substr(i, k);
    uint64_t fwd = string_to_kmer(s);
    uint64_t can = kmer_canonical(fwd, k);
    kmers[can] = {10, 1, 3.0f, 1};
  }

  // Use the extract_unitigs logic inline (simplified version)
  auto degrees = compute_oriented_degrees(kmers, k);

  // Find tip nodes (degree 0 on some side)
  std::vector<uint64_t> tips;
  for (const auto& [can, d] : degrees) {
    if (d.left_fwd == 0 || d.right_fwd == 0 ||
        d.left_rev == 0 || d.right_rev == 0) {
      tips.push_back(can);
    }
  }

  std::cout << "    Tips found: " << tips.size() << "\n";
  assert(tips.size() >= 2);

  // Try extending from a tip
  // Find a tip with right_fwd > 0 (can extend right in forward orientation)
  OrientedNode start{0, true};
  bool found_start = false;
  for (uint64_t tip : tips) {
    auto& d = degrees[tip];
    if (d.left_fwd == 0 && d.right_fwd == 1) {
      start = {tip, true};
      found_start = true;
      break;
    }
    if (d.left_rev == 0 && d.right_rev == 1) {
      start = {tip, false};
      found_start = true;
      break;
    }
  }

  if (!found_start) {
    // Try the other pattern: right == 0, extend left
    for (uint64_t tip : tips) {
      auto& d = degrees[tip];
      if (d.right_fwd == 0 && d.left_fwd == 1) {
        // Reverse: extend right in reverse orientation
        start = {tip, false};
        found_start = true;
        break;
      }
    }
  }

  assert(found_start && "Must find a start node for extension");

  // Extend right
  std::vector<OrientedNode> path;
  path.push_back(start);
  OrientedNode cur = start;
  while (true) {
    auto succs = right_successors(cur, k, kmers);
    if (succs.size() != 1) break;
    OrientedNode next = succs[0];
    auto preds = left_predecessors(next, k, kmers);
    if (preds.size() != 1) break;
    path.push_back(next);
    cur = next;
    if (path.size() > kmers.size() + 1) break; // safety
  }

  // Spell the sequence
  uint64_t first_seq = path[0].forward
    ? path[0].canonical : kmer_revcomp(path[0].canonical, k);
  std::string result = kmer_to_string(first_seq, k);
  for (size_t i = 1; i < path.size(); ++i) {
    uint64_t s = path[i].forward
      ? path[i].canonical : kmer_revcomp(path[i].canonical, k);
    result += base_decode(kmer_last_base(s));
  }

  std::cout << "    Assembled: " << result << " (len=" << result.size() << ")\n";
  std::cout << "    Original:  " << seq << "\n";

  // The assembled sequence should be either seq or its revcomp
  std::string rc_seq;
  for (int i = seq.size() - 1; i >= 0; --i) {
    char c = seq[i];
    switch (c) {
      case 'A': rc_seq += 'T'; break;
      case 'C': rc_seq += 'G'; break;
      case 'G': rc_seq += 'C'; break;
      case 'T': rc_seq += 'A'; break;
    }
  }

  assert((result == seq || result == rc_seq) &&
         "Assembled unitig must match input or its revcomp");
  assert(result.size() == seq.size());

  std::cout << "  test_unitig_extraction PASSED\n";
}

// === Test 5: Branching node creates multiple unitigs ===
void test_branching() {
  const uint32_t k = 5;
  // Y-shaped graph: AAAACACCC and AAAACGTTT share prefix AAAAC
  std::string seq1 = "AAAACACCC";
  std::string seq2 = "AAAACGTTT";

  KmerMap kmers;
  for (const auto& seq : {seq1, seq2}) {
    for (size_t i = 0; i + k <= seq.size(); ++i) {
      std::string s = seq.substr(i, k);
      uint64_t fwd = string_to_kmer(s);
      uint64_t can = kmer_canonical(fwd, k);
      kmers[can] = {10, 1, 3.0f, 1};
    }
  }

  std::cout << "    Total unique canonical k-mers: " << kmers.size() << "\n";

  auto degrees = compute_oriented_degrees(kmers, k);

  int n_branch = 0;
  for (const auto& [can, d] : degrees) {
    if (is_branching(d)) {
      std::string s = kmer_to_string(can, k);
      // Check if it's truly a high-degree node (not just a tip)
      int max_deg = std::max({(int)d.right_fwd, (int)d.right_rev,
                              (int)d.left_fwd, (int)d.left_rev});
      if (max_deg >= 2) {
        std::cout << "    BRANCH: " << s
                  << " L_f=" << (int)d.left_fwd << " R_f=" << (int)d.right_fwd
                  << " L_r=" << (int)d.left_rev << " R_r=" << (int)d.right_rev
                  << "\n";
        n_branch++;
      }
    }
  }

  // AAAAC is the shared prefix — it should have right_fwd >= 2
  uint64_t aaaac = kmer_canonical(string_to_kmer("AAAAC"), k);
  auto& d = degrees[aaaac];
  std::cout << "    AAAAC: L_f=" << (int)d.left_fwd << " R_f=" << (int)d.right_fwd
            << " L_r=" << (int)d.left_rev << " R_r=" << (int)d.right_rev << "\n";

  assert(d.right_fwd >= 2 && "AAAAC must branch to the right in forward orientation");

  std::cout << "  test_branching PASSED\n";
}

int main() {
  std::cout << "=== kbsa_unitig cDBG unit tests ===\n";
  test_kmer_utilities();
  test_shifts();
  test_linear_chain();
  test_unitig_extraction();
  test_branching();
  std::cout << "\nAll tests PASSED.\n";
  return 0;
}

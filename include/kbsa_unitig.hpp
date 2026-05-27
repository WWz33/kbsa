#pragma once

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>

#include "kbsa_core.hpp"
#include "kbsa_params.hpp"
#include "kbsa_merge.hpp"
#include "kbsa_anchor.hpp"

namespace kbsa {

struct KmerInfo {
  uint32_t count1;
  uint32_t count2;
  float z_score;
  int8_t direction;  // +1 = BULK1-enriched, -1 = BULK2-enriched
};

using KmerMap = std::unordered_map<uint64_t, KmerInfo>;

// --- Canonical k-mer utilities (k <= 31, 2-bit packed into uint64_t) ---

inline uint64_t kmer_revcomp(uint64_t kmer, uint32_t k) noexcept
{
  uint64_t rc = 0;
  for (uint32_t i = 0; i < k; ++i) {
    rc = (rc << 2) | (3 - (kmer & 3));
    kmer >>= 2;
  }
  return rc;
}

inline uint64_t kmer_canonical(uint64_t kmer, uint32_t k) noexcept
{
  uint64_t rc = kmer_revcomp(kmer, k);
  return std::min(kmer, rc);
}

inline uint64_t kmer_shift_right(uint64_t kmer, uint32_t k, uint8_t base) noexcept
{
  uint64_t mask = (k < 32) ? ((uint64_t(1) << (2 * k)) - 1) : ~uint64_t(0);
  return ((kmer << 2) | base) & mask;
}

inline uint64_t kmer_shift_left(uint64_t kmer, uint32_t k, uint8_t base) noexcept
{
  return (kmer >> 2) | (uint64_t(base) << (2 * (k - 1)));
}

inline uint8_t kmer_last_base(uint64_t kmer) noexcept
{
  return kmer & 3;
}

inline uint8_t kmer_first_base(uint64_t kmer, uint32_t k) noexcept
{
  return (kmer >> (2 * (k - 1))) & 3;
}

// Encode DNA char to 2-bit (A=0, C=1, G=2, T=3). Non-ACGT → 0.
inline uint8_t base_encode(char c) noexcept
{
  static const auto kBase2 = []() {
    std::array<uint8_t, 256> t {};
    t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
    t['a'] = 0; t['c'] = 1; t['g'] = 2; t['t'] = 3;
    return t;
  }();
  return kBase2[static_cast<unsigned char>(c)];
}

inline char base_decode(uint8_t b) noexcept
{
  return "ACGT"[b & 3];
}

inline std::string kmer_to_string(uint64_t kmer, uint32_t k)
{
  std::string s(k, 'A');
  for (int i = k - 1; i >= 0; --i) {
    s[i] = base_decode(kmer & 3);
    kmer >>= 2;
  }
  return s;
}

inline uint64_t string_to_kmer(const std::string& s)
{
  uint64_t kmer = 0;
  for (char c : s) {
    kmer = (kmer << 2) | base_encode(c);
  }
  return kmer;
}

// --- Oriented neighbor enumeration ---

struct OrientedNode {
  uint64_t canonical;
  bool forward;  // true = canonical == sequence as read
};

// Given an oriented node, enumerate right successors that exist in the map.
// Returns oriented nodes for each valid successor.
inline std::vector<OrientedNode> right_successors(
    OrientedNode node, uint32_t k, const KmerMap& kmers)
{
  // Recover the actual sequence k-mer
  uint64_t seq_kmer = node.forward
    ? node.canonical
    : kmer_revcomp(node.canonical, k);

  std::vector<OrientedNode> result;
  for (uint8_t b = 0; b < 4; ++b) {
    uint64_t candidate = kmer_shift_right(seq_kmer, k, b);
    uint64_t can = kmer_canonical(candidate, k);
    if (kmers.count(can)) {
      bool fwd = (candidate == can);
      result.push_back({can, fwd});
    }
  }
  return result;
}

// Given an oriented node, enumerate left predecessors that exist in the map.
inline std::vector<OrientedNode> left_predecessors(
    OrientedNode node, uint32_t k, const KmerMap& kmers)
{
  uint64_t seq_kmer = node.forward
    ? node.canonical
    : kmer_revcomp(node.canonical, k);

  std::vector<OrientedNode> result;
  for (uint8_t b = 0; b < 4; ++b) {
    uint64_t candidate = kmer_shift_left(seq_kmer, k, b);
    uint64_t can = kmer_canonical(candidate, k);
    if (kmers.count(can)) {
      bool fwd = (candidate == can);
      result.push_back({can, fwd});
    }
  }
  return result;
}

// --- Degree computation ---

struct NodeDegree {
  uint8_t left_fwd;
  uint8_t right_fwd;
  uint8_t left_rev;
  uint8_t right_rev;
};

inline std::unordered_map<uint64_t, NodeDegree> compute_oriented_degrees(
    const KmerMap& kmers, uint32_t k)
{
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

inline bool is_branching(const NodeDegree& d) noexcept
{
  // Branching if not simple pass-through in either orientation
  bool simple_fwd = (d.left_fwd == 1 && d.right_fwd == 1);
  bool simple_rev = (d.left_rev == 1 && d.right_rev == 1);
  return !simple_fwd || !simple_rev;
}

// --- Differential k-mer extraction ---

struct DiffThresholds {
  uint32_t enriched_min;
  uint32_t opposite_max;
};

inline uint64_t find_valley_depth(const std::vector<HistEntry>& hist)
{
  if (hist.empty()) return 3;
  // hist is freq-ascending (KMC kmc_tools histogram order).
  size_t start = 0;
  for (size_t i = 0; i < hist.size(); ++i) {
    if (hist[i].freq >= 2) { start = i; break; }
  }
  for (size_t i = start + 1; i < hist.size(); ++i) {
    if (hist[i].count > hist[i - 1].count)
      return hist[i - 1].freq;
  }
  return 3;
}

inline DiffThresholds derive_thresholds(const std::vector<HistEntry>& hist)
{
  uint64_t peak = find_peak_depth(hist);
  uint64_t valley = find_valley_depth(hist);
  // kcov ≈ peak/2 for diploid (peak is homozygous = 2×haploid coverage)
  uint64_t kcov = std::max<uint64_t>(peak / 2, 1);
  uint64_t kcov_guard = std::max<uint64_t>(3, (kcov + 2) / 3); // ceil(kcov/3)

  DiffThresholds t {};
  t.enriched_min = static_cast<uint32_t>(std::max(valley, kcov_guard));
  t.opposite_max = static_cast<uint32_t>(
      std::max<uint64_t>(1, std::min(valley / 2, (kcov + 5) / 6))); // min(valley/2, ceil(kcov/6))
  return t;
}

// Convert string k-mer to packed uint64 canonical form
inline uint64_t string_kmer_to_canonical(const char* s, uint32_t k)
{
  uint64_t fwd = 0;
  for (uint32_t i = 0; i < k; ++i)
    fwd = (fwd << 2) | base_encode(s[i]);
  return kmer_canonical(fwd, k);
}

inline KmerMap extract_differential_kmers(
    const std::string& bulk1_db_path,
    const std::string& bulk2_db_path,
    const std::vector<HistEntry>& hist1,
    const std::vector<HistEntry>& hist2,
    const Params& params,
    uint32_t enriched_min_override = 0,
    uint32_t opposite_max_override = 0)
{
  DiffThresholds t1 = derive_thresholds(hist1);
  DiffThresholds t2 = derive_thresholds(hist2);

  uint32_t enr_min1 = enriched_min_override ? enriched_min_override : t1.enriched_min;
  uint32_t enr_min2 = enriched_min_override ? enriched_min_override : t2.enriched_min;
  uint32_t opp_max1 = opposite_max_override ? opposite_max_override : t1.opposite_max;
  uint32_t opp_max2 = opposite_max_override ? opposite_max_override : t2.opposite_max;

  KmerMap result;
  KmcMergeIterator iter(bulk1_db_path, bulk2_db_path);
  uint32_t k = iter.kmer_length();

  if (k > 31)
    throw std::runtime_error("kbsa unitig requires k <= 31, got k=" + std::to_string(k));

  MergedKmer mk;
  while (iter.next(mk)) {
    uint32_t c1 = static_cast<uint32_t>(mk.bulk1_count);
    uint32_t c2 = static_cast<uint32_t>(mk.bulk2_count);

    bool bulk1_enriched = (c1 >= enr_min1 && c2 <= opp_max2);
    bool bulk2_enriched = (c2 >= enr_min2 && c1 <= opp_max1);

    if (!bulk1_enriched && !bulk2_enriched) continue;

    Score sc = score_kmer(c2, c1,
                          params.bulk1_to_bulk2_scale,
                          0, UINT64_MAX,
                          0.0, 0.0,
                          params.error_rate);

    KmerInfo info {};
    info.count1 = c1;
    info.count2 = c2;
    info.z_score = static_cast<float>(sc.z_score);
    info.direction = (sc.sign == Significance::BULK1) ? int8_t(+1) : int8_t(-1);

    uint64_t can = string_kmer_to_canonical(mk.kmer_str, k);
    result[can] = info;
  }

  return result;
}

inline KmerMap extract_differential_kmers(
    const std::string& bulk1_db_path,
    const std::string& bulk2_db_path,
    const std::string& bulk1_hist_path,
    const std::string& bulk2_hist_path,
    const Params& params,
    uint32_t enriched_min_override = 0,
    uint32_t opposite_max_override = 0)
{
  auto h1 = load_histogram(bulk1_hist_path);
  auto h2 = load_histogram(bulk2_hist_path);
  return extract_differential_kmers(
      bulk1_db_path, bulk2_db_path, h1, h2,
      params, enriched_min_override, opposite_max_override);
}

// --- Unitig extraction (edge-based traversal) ---

inline std::vector<std::string> extract_unitigs(
    const KmerMap& kmers, uint32_t k, uint32_t min_length)
{
  if (kmers.empty()) return {};

  auto degrees = compute_oriented_degrees(kmers, k);

  // Track visited directed edges as (from_canonical, to_canonical) pairs.
  // After collecting unitigs, deduplicate by removing reverse complements.
  struct PairHash {
    size_t operator()(std::pair<uint64_t,uint64_t> p) const noexcept {
      return std::hash<uint64_t>{}(p.first) ^ (std::hash<uint64_t>{}(p.second) << 32);
    }
  };
  std::unordered_set<std::pair<uint64_t,uint64_t>, PairHash> visited_edges;

  auto make_edge = [](uint64_t a, uint64_t b) {
    return std::make_pair(a, b);
  };

  std::vector<std::string> unitigs;

  // Find starting points: nodes that are branching or tips
  std::vector<uint64_t> starts;
  for (const auto& [can, deg] : degrees) {
    if (is_branching(deg) || deg.left_fwd == 0 || deg.right_fwd == 0 ||
        deg.left_rev == 0 || deg.right_rev == 0) {
      starts.push_back(can);
    }
  }

  auto extend_right = [&](OrientedNode start) -> std::vector<OrientedNode> {
    std::vector<OrientedNode> path;
    path.push_back(start);
    OrientedNode cur = start;

    while (true) {
      auto succs = right_successors(cur, k, kmers);
      if (succs.size() != 1) break;

      OrientedNode next = succs[0];
      auto edge = make_edge(cur.canonical, next.canonical);
      if (visited_edges.count(edge)) break;

      // Check if next is simple pass-through from this direction
      auto preds = left_predecessors(next, k, kmers);
      if (preds.size() != 1) break;

      // Mark both directions (same biological edge in bidirected graph)
      visited_edges.insert(edge);
      visited_edges.insert(make_edge(next.canonical, cur.canonical));
      path.push_back(next);
      cur = next;
    }
    return path;
  };

  auto path_to_sequence = [&](const std::vector<OrientedNode>& path) -> std::string {
    if (path.empty()) return "";
    // First node contributes full k-mer
    uint64_t first_seq = path[0].forward
      ? path[0].canonical
      : kmer_revcomp(path[0].canonical, k);
    std::string seq = kmer_to_string(first_seq, k);
    // Each subsequent node contributes 1 base
    for (size_t i = 1; i < path.size(); ++i) {
      uint64_t s = path[i].forward
        ? path[i].canonical
        : kmer_revcomp(path[i].canonical, k);
      seq += base_decode(kmer_last_base(s));
    }
    return seq;
  };

  // Traverse from each starting node in both orientations
  for (uint64_t start_can : starts) {
    for (bool fwd : {true, false}) {
      OrientedNode start {start_can, fwd};
      auto succs = right_successors(start, k, kmers);

      for (auto& first_succ : succs) {
        auto edge = make_edge(start_can, first_succ.canonical);
        if (visited_edges.count(edge)) continue;
        visited_edges.insert(edge);
        visited_edges.insert(make_edge(first_succ.canonical, start_can));

        // Extend right from first_succ
        auto right_path = extend_right(first_succ);

        // Full path = start + right_path
        std::vector<OrientedNode> full_path;
        full_path.push_back(start);
        full_path.insert(full_path.end(), right_path.begin(), right_path.end());

        std::string seq = path_to_sequence(full_path);
        if (seq.size() >= min_length) {
          unitigs.push_back(std::move(seq));
        }
      }
    }
  }

  // Handle circular components: find any unvisited k-mers
  std::unordered_set<uint64_t> in_unitig;
  for (const auto& [can, deg] : degrees) {
    if (is_branching(deg) || deg.left_fwd == 0 || deg.right_fwd == 0 ||
        deg.left_rev == 0 || deg.right_rev == 0) {
      in_unitig.insert(can);
    }
  }
  // Nodes already in paths via edges
  for (const auto& [from, to] : visited_edges) {
    in_unitig.insert(from);
    in_unitig.insert(to);
  }

  for (const auto& [can, _] : kmers) {
    if (in_unitig.count(can)) continue;
    // This node is in a cycle — break it here
    OrientedNode start {can, true};
    auto right_path = extend_right(start);
    if (right_path.size() > 1) {
      std::string seq = path_to_sequence(right_path);
      if (seq.size() >= min_length) {
        unitigs.push_back(std::move(seq));
      }
    }
    for (auto& n : right_path) in_unitig.insert(n.canonical);
  }

  // Deduplicate: keep canonical form (lexicographically smaller of seq and revcomp)
  auto seq_revcomp = [](const std::string& s) -> std::string {
    std::string rc(s.size(), 'A');
    for (size_t i = 0; i < s.size(); ++i) {
      char c = s[s.size() - 1 - i];
      switch (c) {
        case 'A': rc[i] = 'T'; break;
        case 'C': rc[i] = 'G'; break;
        case 'G': rc[i] = 'C'; break;
        case 'T': rc[i] = 'A'; break;
        default:  rc[i] = 'N'; break;
      }
    }
    return rc;
  };

  std::unordered_set<std::string> seen;
  std::vector<std::string> deduped;
  for (auto& u : unitigs) {
    std::string rc = seq_revcomp(u);
    std::string canonical_seq = std::min(u, rc);
    if (!seen.count(canonical_seq)) {
      seen.insert(canonical_seq);
      deduped.push_back(std::move(u));
    }
  }

  return deduped;
}

// --- Two-phase extension: extend unitigs using full KMC databases ---

inline std::string extend_unitig_right(
    const std::string& unitig, uint32_t k,
    KmcRandomAccess& db1, KmcRandomAccess& db2,
    uint32_t min_ext_count = 2, uint32_t max_ext = 10000)
{
  std::string seq = unitig;
  std::unordered_set<uint64_t> visited;
  char buf[33];
  buf[k] = '\0';

  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = 0;
    for (uint32_t j = 0; j < k; ++j)
      fwd = (fwd << 2) | base_encode(seq[i + j]);
    visited.insert(kmer_canonical(fwd, k));
  }

  for (uint32_t ext = 0; ext < max_ext; ++ext) {
    std::memcpy(buf, seq.data() + seq.size() - k + 1, k - 1);
    int n_succ = 0;
    char succ_base = 'N';

    for (char b : {'A', 'C', 'G', 'T'}) {
      buf[k - 1] = b;
      uint64_t fwd = 0;
      for (uint32_t j = 0; j < k; ++j)
        fwd = (fwd << 2) | base_encode(buf[j]);
      uint64_t can = kmer_canonical(fwd, k);
      if (visited.count(can)) continue;

      uint64_t c1 = db1.query(buf);
      uint64_t c2 = db2.query(buf);
      if (c1 + c2 >= min_ext_count) {
        ++n_succ;
        succ_base = b;
      }
    }

    if (n_succ != 1) break;

    seq += succ_base;
    uint64_t new_fwd = 0;
    for (uint32_t j = 0; j < k; ++j)
      new_fwd = (new_fwd << 2) | base_encode(seq[seq.size() - k + j]);
    visited.insert(kmer_canonical(new_fwd, k));
  }
  return seq;
}

inline std::string extend_unitig_left(
    const std::string& unitig, uint32_t k,
    KmcRandomAccess& db1, KmcRandomAccess& db2,
    uint32_t min_ext_count = 2, uint32_t max_ext = 10000)
{
  std::string seq = unitig;
  std::unordered_set<uint64_t> visited;
  char buf[33];
  buf[k] = '\0';

  for (size_t i = 0; i + k <= seq.size(); ++i) {
    uint64_t fwd = 0;
    for (uint32_t j = 0; j < k; ++j)
      fwd = (fwd << 2) | base_encode(seq[i + j]);
    visited.insert(kmer_canonical(fwd, k));
  }

  for (uint32_t ext = 0; ext < max_ext; ++ext) {
    std::memcpy(buf + 1, seq.data(), k - 1);
    int n_pred = 0;
    char pred_base = 'N';

    for (char b : {'A', 'C', 'G', 'T'}) {
      buf[0] = b;
      uint64_t fwd = 0;
      for (uint32_t j = 0; j < k; ++j)
        fwd = (fwd << 2) | base_encode(buf[j]);
      uint64_t can = kmer_canonical(fwd, k);
      if (visited.count(can)) continue;

      uint64_t c1 = db1.query(buf);
      uint64_t c2 = db2.query(buf);
      if (c1 + c2 >= min_ext_count) {
        ++n_pred;
        pred_base = b;
      }
    }

    if (n_pred != 1) break;

    seq = pred_base + seq;
    uint64_t new_fwd = 0;
    for (uint32_t j = 0; j < k; ++j)
      new_fwd = (new_fwd << 2) | base_encode(seq[j]);
    visited.insert(kmer_canonical(new_fwd, k));
  }
  return seq;
}

inline std::vector<std::string> extend_unitigs_two_phase(
    const std::vector<std::string>& unitigs, uint32_t k,
    const std::string& bulk1_db_path,
    const std::string& bulk2_db_path,
    uint32_t min_ext_count = 2)
{
  KmcRandomAccess db1, db2;
  if (!db1.open(bulk1_db_path))
    throw std::runtime_error("Cannot open KMC RA: " + bulk1_db_path);
  if (!db2.open(bulk2_db_path))
    throw std::runtime_error("Cannot open KMC RA: " + bulk2_db_path);

  std::vector<std::string> extended;
  extended.reserve(unitigs.size());

  for (const auto& u : unitigs) {
    std::string ext = extend_unitig_right(u, k, db1, db2, min_ext_count);
    ext = extend_unitig_left(ext, k, db1, db2, min_ext_count);
    extended.push_back(std::move(ext));
  }

  db1.close();
  db2.close();
  return extended;
}

// --- BCALM2 pipeline: stream differential k-mers to FASTA, no in-memory KmerMap ---

struct DiffStats {
  uint64_t total_diff_kmers;
  uint32_t k;
  uint32_t enriched_min;
  uint32_t opposite_max;
};

inline DiffStats stream_differential_kmers_to_fasta(
    const std::string& bulk1_db_path,
    const std::string& bulk2_db_path,
    const std::vector<HistEntry>& hist1,
    const std::vector<HistEntry>& hist2,
    const std::string& fasta_out_path,
    uint32_t enriched_min_override = 0,
    uint32_t opposite_max_override = 0,
    const std::string& parent1_db_path = "",
    const std::string& parent2_db_path = "",
    uint32_t parent_min_count = 3)
{
  DiffThresholds t1 = derive_thresholds(hist1);
  DiffThresholds t2 = derive_thresholds(hist2);
  uint32_t enr_min1 = enriched_min_override ? enriched_min_override : t1.enriched_min;
  uint32_t enr_min2 = enriched_min_override ? enriched_min_override : t2.enriched_min;
  uint32_t opp_max1 = opposite_max_override ? opposite_max_override : t1.opposite_max;
  uint32_t opp_max2 = opposite_max_override ? opposite_max_override : t2.opposite_max;

  // Optional parental segregating filter
  bool use_parents = !parent1_db_path.empty() && !parent2_db_path.empty();
  KmcRandomAccess p1_db, p2_db;
  if (use_parents) {
    if (!p1_db.open(parent1_db_path))
      throw std::runtime_error("Cannot open parent1 KMC RA: " + parent1_db_path);
    if (!p2_db.open(parent2_db_path))
      throw std::runtime_error("Cannot open parent2 KMC RA: " + parent2_db_path);
  }

  KmcMergeIterator iter(bulk1_db_path, bulk2_db_path);
  uint32_t k = iter.kmer_length();
  if (k > 31)
    throw std::runtime_error("kbsa unitig requires k <= 31, got k=" + std::to_string(k));

  std::ofstream fa(fasta_out_path);
  if (!fa)
    throw std::runtime_error("Cannot open output FASTA: " + fasta_out_path);
  static char buf[1 << 20];
  fa.rdbuf()->pubsetbuf(buf, sizeof(buf));

  uint64_t total = 0, filtered_by_parents = 0;
  MergedKmer mk;
  while (iter.next(mk)) {
    uint32_t c1 = static_cast<uint32_t>(mk.bulk1_count);
    uint32_t c2 = static_cast<uint32_t>(mk.bulk2_count);
    bool bulk1_enriched = (c1 >= enr_min1 && c2 <= opp_max2);
    bool bulk2_enriched = (c2 >= enr_min2 && c1 <= opp_max1);
    if (!bulk1_enriched && !bulk2_enriched) continue;

    if (use_parents) {
      uint64_t p1c = p1_db.query(mk.kmer_str);
      uint64_t p2c = p2_db.query(mk.kmer_str);
      bool p1_present = (p1c >= parent_min_count);
      bool p2_present = (p2c >= parent_min_count);
      bool segregating = (p1_present != p2_present);
      if (!segregating) { ++filtered_by_parents; continue; }
    }

    fa << ">k" << total << '\n' << mk.kmer_str << '\n';
    ++total;
  }
  fa.close();
  if (use_parents) {
    p1_db.close();
    p2_db.close();
    fprintf(stderr, "[unitig] Parental filter: kept %lu, removed %lu non-segregating k-mers\n",
            total, filtered_by_parents);
  }

  DiffStats s {};
  s.total_diff_kmers = total;
  s.k = k;
  s.enriched_min = std::max(enr_min1, enr_min2);
  s.opposite_max = std::max(opp_max1, opp_max2);
  return s;
}

// Read unitigs from BCALM2 output FASTA
inline std::vector<std::string> read_bcalm2_unitigs(
    const std::string& fasta_path, uint32_t min_length)
{
  std::vector<std::string> unitigs;
  std::ifstream fa(fasta_path);
  if (!fa) throw std::runtime_error("Cannot open BCALM2 output: " + fasta_path);

  std::string line, seq;
  while (std::getline(fa, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (!seq.empty() && seq.size() >= min_length)
        unitigs.push_back(std::move(seq));
      seq.clear();
    } else {
      seq += line;
    }
  }
  if (!seq.empty() && seq.size() >= min_length)
    unitigs.push_back(std::move(seq));
  return unitigs;
}

}  // namespace kbsa

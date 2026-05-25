#pragma once

#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include "kbsa_unitig.hpp"

namespace kbsa {

struct MinimizerHit {
  uint64_t ref_pos;
  uint32_t unitig_pos;
  bool same_strand;
};

struct MappingResult {
  std::string chrom;
  uint64_t start;
  uint64_t end;
  bool forward;
  uint32_t n_anchors;
  double score;
};

inline uint64_t minimizer_hash(const char* seq, uint32_t len) noexcept
{
  uint64_t fwd = 0, rev = 0;
  for (uint32_t i = 0; i < len; ++i) {
    uint8_t b = base_encode(seq[i]);
    fwd = (fwd << 2) | b;
    rev = (rev >> 2) | (uint64_t(3 - b) << (2 * (len - 1)));
  }
  return std::min(fwd, rev);
}

inline double seq_entropy(const char* seq, uint32_t len) noexcept
{
  uint32_t counts[4] = {};
  for (uint32_t i = 0; i < len; ++i)
    counts[base_encode(seq[i])]++;
  double ent = 0.0;
  for (int i = 0; i < 4; ++i) {
    if (counts[i] == 0) continue;
    double p = static_cast<double>(counts[i]) / len;
    ent -= p * std::log2(p);
  }
  return ent;
}

struct RefMinimizerIndex {
  struct RefPos {
    uint32_t chrom_idx;
    uint64_t pos;
  };

  std::unordered_map<uint64_t, std::vector<RefPos>> index;
  std::vector<std::string> chrom_names;
  uint32_t mini_k;
  uint32_t freq_cap;
  double min_entropy;

  RefMinimizerIndex(uint32_t k = 15, uint32_t cap = 50, double ent = 1.5)
    : mini_k(k), freq_cap(cap), min_entropy(ent) {}

  void build(const std::vector<std::pair<std::string, std::string>>& chroms)
  {
    chrom_names.clear();
    index.clear();

    for (size_t ci = 0; ci < chroms.size(); ++ci) {
      chrom_names.push_back(chroms[ci].first);
      const std::string& seq = chroms[ci].second;
      if (seq.size() < mini_k) continue;

      for (size_t i = 0; i + mini_k <= seq.size(); ++i) {
        const char* s = seq.data() + i;

        bool has_n = false;
        for (uint32_t j = 0; j < mini_k; ++j) {
          if (s[j] == 'N' || s[j] == 'n') { has_n = true; break; }
        }
        if (has_n) continue;

        if (seq_entropy(s, mini_k) < min_entropy) continue;

        uint64_t h = minimizer_hash(s, mini_k);
        index[h].push_back({static_cast<uint32_t>(ci), i});
      }
    }

    for (auto it = index.begin(); it != index.end(); ) {
      if (it->second.size() > freq_cap)
        it = index.erase(it);
      else
        ++it;
    }
  }
};


// LIS-based co-linear chaining — returns the actual LIS subsequence
inline std::vector<uint64_t> lis_chain(const std::vector<uint64_t>& positions)
{
  if (positions.empty()) return {};
  size_t n = positions.size();
  std::vector<uint64_t> tails;
  std::vector<size_t> tail_idx;
  std::vector<size_t> parent(n, SIZE_MAX);
  std::vector<size_t> idx_in_tails(n);

  for (size_t i = 0; i < n; ++i) {
    auto it = std::lower_bound(tails.begin(), tails.end(), positions[i]);
    size_t pos_in_tails = static_cast<size_t>(it - tails.begin());
    if (it == tails.end()) {
      tails.push_back(positions[i]);
      tail_idx.push_back(i);
    } else {
      *it = positions[i];
      tail_idx[pos_in_tails] = i;
    }
    idx_in_tails[i] = pos_in_tails;
    if (pos_in_tails > 0) parent[i] = tail_idx[pos_in_tails - 1];
  }

  std::vector<uint64_t> result(tails.size());
  size_t cur = tail_idx.back();
  for (int j = static_cast<int>(tails.size()) - 1; j >= 0; --j) {
    result[j] = positions[cur];
    cur = parent[cur];
  }
  return result;
}

inline MappingResult map_unitig(
    const std::string& unitig, const RefMinimizerIndex& idx,
    uint32_t min_anchors = 3)
{
  MappingResult best {};
  best.n_anchors = 0;
  uint64_t max_span = unitig.size() * 2;

  // Collect all hits grouped by chromosome
  std::unordered_map<uint32_t, std::vector<MinimizerHit>> hits_by_chrom;
  if (unitig.size() >= idx.mini_k) {
    for (size_t i = 0; i + idx.mini_k <= unitig.size(); ++i) {
      const char* s = unitig.data() + i;
      if (seq_entropy(s, idx.mini_k) < idx.min_entropy) continue;
      uint64_t h = minimizer_hash(s, idx.mini_k);
      auto it = idx.index.find(h);
      if (it == idx.index.end()) continue;
      for (const auto& rp : it->second) {
        hits_by_chrom[rp.chrom_idx].push_back(
            {rp.pos, static_cast<uint32_t>(i), true});
      }
    }
  }

  for (auto& [ci, hits] : hits_by_chrom) {
    if (hits.size() < min_anchors) continue;

    std::sort(hits.begin(), hits.end(),
      [](const MinimizerHit& a, const MinimizerHit& b) {
        return a.unitig_pos < b.unitig_pos;
      });

    std::vector<uint64_t> fwd_pos;
    for (const auto& h : hits) fwd_pos.push_back(h.ref_pos);

    auto chain = lis_chain(fwd_pos);
    if (chain.size() < min_anchors) continue;

    uint32_t best_window_len = 0;
    uint64_t best_window_start = 0, best_window_end = 0;
    size_t left = 0;
    for (size_t right = 0; right < chain.size(); ++right) {
      while (chain[right] - chain[left] > max_span) ++left;
      uint32_t wlen = static_cast<uint32_t>(right - left + 1);
      if (wlen > best_window_len) {
        best_window_len = wlen;
        best_window_start = chain[left];
        best_window_end = chain[right] + idx.mini_k;
      }
    }

    if (best_window_len >= min_anchors && best_window_len > best.n_anchors) {
      best.chrom = idx.chrom_names[ci];
      best.n_anchors = best_window_len;
      best.forward = true;
      best.start = best_window_start;
      best.end = best_window_end;
    }
  }

  return best;
}

}  // namespace kbsa

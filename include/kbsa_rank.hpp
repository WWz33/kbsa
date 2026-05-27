#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace kbsa {

// Pack canonical 31-mer into uint64 (2 bits per base, 62 bits used).
inline uint64_t pack31(const char* s) noexcept
{
  // 0xFF = invalid base; A=0, C=1, G=2, T=3.
  static const auto kBase2 = []() {
    std::array<uint8_t, 256> t {};
    t.fill(0xFF);
    t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
    return t;
  }();
  uint64_t x = 0;
  for (int i = 0; i < 31; ++i) {
    uint8_t c = kBase2[static_cast<unsigned char>(s[i])];
    if (c == 0xFF) return UINT64_MAX;
    x = (x << 2) | c;
  }
  return x;
}

inline uint64_t reverse_complement_packed(uint64_t x, int k = 31) noexcept
{
  uint64_t rc = 0;
  for (int i = 0; i < k; ++i) {
    uint64_t base = x & 3ULL;
    rc = (rc << 2) | (base ^ 3ULL);
    x >>= 2;
  }
  return rc;
}

inline uint64_t canonical_packed(uint64_t fwd) noexcept
{
  uint64_t rc = reverse_complement_packed(fwd, 31);
  return fwd < rc ? fwd : rc;
}

// Per-k-mer score loaded from kbsa score TSV.
struct KmerEntry {
  float z_score;     // signed: + CONTROL-enriched, - CASE-enriched
  float kai_reg;
  uint32_t bulk1_raw;
  uint32_t bulk2_raw;
};

// Per-unitig aggregate.
struct UnitigStats {
  std::string id;
  uint32_t length = 0;
  uint32_t n_kmers_total = 0;     // n_kmers in unitig sequence (length-30)
  uint32_t n_kmers_scored = 0;    // n_kmers that were in score table
  uint32_t n_bulk1 = 0;            // n_kmers with z<0 (CASE-enriched)
  uint32_t n_bulk2 = 0;            // n_kmers with z>0
  double sum_abs_z = 0.0;
  double sum_signed_z = 0.0;
  double max_abs_z = 0.0;
  double median_abs_z = 0.0;
  std::string sequence;
};

// Load kbsa score TSV into a hash table keyed by canonical packed 31-mer.
// Skips header. Filters by min_g (g_score column) for memory efficiency.
inline size_t load_score_table(const std::string& tsv_path,
                               double min_g,
                               std::unordered_map<uint64_t, KmerEntry>& out)
{
  std::ifstream in(tsv_path);
  if (!in.good()) {
    fprintf(stderr, "ERROR: cannot open score TSV: %s\n", tsv_path.c_str());
    return 0;
  }
  std::string line;
  std::getline(in, line);  // header

  size_t n_loaded = 0, n_skipped = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    // columns: kmer  significance  kai_reg  g_score  z_score  bulk2_raw  bulk1_raw  bulk2_adj  bulk1_adj
    size_t p0 = 0;
    auto next_field = [&](size_t& p) -> std::string {
      size_t e = line.find('\t', p);
      std::string s = line.substr(p, e - p);
      p = (e == std::string::npos) ? line.size() : e + 1;
      return s;
    };
    std::string kmer = next_field(p0);
    next_field(p0);  // sig
    double kai = std::stod(next_field(p0));
    double g   = std::stod(next_field(p0));
    if (g < min_g) { ++n_skipped; continue; }
    double z   = std::stod(next_field(p0));
    uint32_t bulk2_raw = static_cast<uint32_t>(std::stoul(next_field(p0)));
    uint32_t bulk1_raw = static_cast<uint32_t>(std::stoul(next_field(p0)));
    if (kmer.size() != 31) continue;
    uint64_t fwd = pack31(kmer.c_str());
    if (fwd == UINT64_MAX) continue;
    uint64_t key = canonical_packed(fwd);
    KmerEntry e;
    e.z_score = static_cast<float>(z);
    e.kai_reg = static_cast<float>(kai);
    e.bulk1_raw = bulk1_raw;
    e.bulk2_raw = bulk2_raw;
    out.emplace(key, e);
    ++n_loaded;
  }
  fprintf(stderr, "[rank] loaded %zu k-mers (skipped %zu below g_min)\n", n_loaded, n_skipped);
  return n_loaded;
}

// Walk a unitig sequence, look up each constituent k-mer, aggregate stats.
inline UnitigStats aggregate_unitig(const std::string& id,
                                    const std::string& seq,
                                    const std::unordered_map<uint64_t, KmerEntry>& tbl)
{
  UnitigStats u;
  u.id = id;
  u.length = static_cast<uint32_t>(seq.size());
  if (seq.size() < 31) return u;
  uint32_t n_total = static_cast<uint32_t>(seq.size() - 31 + 1);
  u.n_kmers_total = n_total;

  std::vector<float> abs_zs;
  abs_zs.reserve(n_total);

  for (size_t i = 0; i + 31 <= seq.size(); ++i) {
    uint64_t fwd = pack31(seq.data() + i);
    if (fwd == UINT64_MAX) continue;
    uint64_t key = canonical_packed(fwd);
    auto it = tbl.find(key);
    if (it == tbl.end()) continue;
    float z = it->second.z_score;
    float az = std::abs(z);
    ++u.n_kmers_scored;
    if (z >= 0) ++u.n_bulk2; else ++u.n_bulk1;
    u.sum_signed_z += z;
    u.sum_abs_z += az;
    if (az > u.max_abs_z) u.max_abs_z = az;
    abs_zs.push_back(az);
  }

  if (!abs_zs.empty()) {
    auto mid = abs_zs.begin() + abs_zs.size() / 2;
    std::nth_element(abs_zs.begin(), mid, abs_zs.end());
    u.median_abs_z = *mid;
  }

  u.sequence = seq;
  return u;
}

}  // namespace kbsa

#pragma once

#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "kbsa_core.hpp"
#include "kbsa_unitig.hpp"

namespace kbsa {

struct UnitigScore {
  std::string id;
  std::string sequence;
  uint32_t length;
  uint32_t n_kmers;
  int8_t dominant_dir;    // +1 = BULK1, -1 = BULK2
  float coherence;        // [0.5, 1.0]
  float mean_z;
  float excess_per_kb;
  float sum_excess_mass;
  uint64_t sum_bulk1;
  uint64_t sum_bulk2;
  float rank_score;       // excess_per_kb * coherence
};

inline UnitigScore score_unitig(const std::string& seq, uint32_t k,
                                const KmerMap& kmers, const Params& params,
                                float tau)
{
  UnitigScore us {};
  us.length = static_cast<uint32_t>(seq.size());
  us.n_kmers = (seq.size() >= k) ? static_cast<uint32_t>(seq.size() - k + 1) : 0;

  if (us.n_kmers == 0) return us;

  int bulk1_count_dir = 0;
  double sum_abs_z = 0.0;
  double sum_excess = 0.0;
  uint64_t total_b1 = 0, total_b2 = 0;

  for (size_t i = 0; i + k <= seq.size(); ++i) {
    // Compute canonical k-mer directly from char* (no substr allocation).
    uint64_t fwd = 0;
    for (uint32_t j = 0; j < k; ++j)
      fwd = (fwd << 2) | base_encode(seq[i + j]);
    uint64_t can = kmer_canonical(fwd, k);

    auto it = kmers.find(can);
    if (it == kmers.end()) continue;

    const KmerInfo& info = it->second;
    total_b1 += info.count1;
    total_b2 += info.count2;

    Score sc = score_kmer(info.count2, info.count1,
                          params.bulk1_to_bulk2_scale,
                          0, UINT64_MAX,
                          0.0, 0.0,
                          params.error_rate);

    double z = sc.z_score;
    sum_abs_z += std::abs(z);

    if (z > 0) bulk1_count_dir++;
    else if (z < 0) bulk1_count_dir--;

    double excess = std::max(0.0, std::abs(z) - tau);
    sum_excess += excess;
  }

  us.sum_bulk1 = total_b1;
  us.sum_bulk2 = total_b2;
  us.dominant_dir = (bulk1_count_dir >= 0) ? int8_t(+1) : int8_t(-1);

  int n_dominant = (us.dominant_dir > 0)
    ? (us.n_kmers + bulk1_count_dir) / 2
    : (us.n_kmers - bulk1_count_dir) / 2;
  us.coherence = static_cast<float>(n_dominant) / us.n_kmers;

  us.mean_z = static_cast<float>(sum_abs_z / us.n_kmers);
  us.sum_excess_mass = static_cast<float>(sum_excess);
  us.excess_per_kb = static_cast<float>(sum_excess / (us.length / 1000.0));
  us.rank_score = us.excess_per_kb * us.coherence;

  return us;
}

inline std::vector<UnitigScore> score_and_rank_unitigs(
    const std::vector<std::string>& unitigs, uint32_t k,
    const KmerMap& kmers, const Params& params, float tau)
{
  std::vector<UnitigScore> scores;
  scores.reserve(unitigs.size());

  for (size_t i = 0; i < unitigs.size(); ++i) {
    auto us = score_unitig(unitigs[i], k, kmers, params, tau);
    us.id = "unitig_" + std::to_string(i + 1);
    us.sequence = unitigs[i];
    scores.push_back(std::move(us));
  }

  std::sort(scores.begin(), scores.end(), [](const UnitigScore& a, const UnitigScore& b) {
    if (a.rank_score != b.rank_score) return a.rank_score > b.rank_score;
    return a.sum_excess_mass > b.sum_excess_mass;
  });

  // Reassign IDs by rank
  for (size_t i = 0; i < scores.size(); ++i) {
    const char* dir_str = (scores[i].dominant_dir > 0) ? "BULK1" : "BULK2";
    scores[i].id = "unitig_" + std::to_string(i + 1) + "_" + dir_str;
  }

  return scores;
}

// --- KMC-based scoring (no KmerMap needed) ---

inline UnitigScore score_unitig_kmc(const std::string& seq, uint32_t k,
                                    KmcRandomAccess& db1, KmcRandomAccess& db2,
                                    const Params& params, float tau)
{
  UnitigScore us {};
  us.length = static_cast<uint32_t>(seq.size());
  us.n_kmers = (seq.size() >= k) ? static_cast<uint32_t>(seq.size() - k + 1) : 0;
  if (us.n_kmers == 0) return us;

  int bulk1_count_dir = 0;
  double sum_abs_z = 0.0, sum_excess = 0.0;
  uint64_t total_b1 = 0, total_b2 = 0;
  char buf[33];
  buf[k] = '\0';

  for (size_t i = 0; i + k <= seq.size(); ++i) {
    std::memcpy(buf, seq.data() + i, k);
    uint64_t c1 = db1.query(buf);
    uint64_t c2 = db2.query(buf);
    total_b1 += c1;
    total_b2 += c2;

    Score sc = score_kmer(c2, c1,
                          params.bulk1_to_bulk2_scale,
                          0, UINT64_MAX, 0.0, 0.0,
                          params.error_rate);
    double z = sc.z_score;
    sum_abs_z += std::abs(z);
    if (z > 0) bulk1_count_dir++;
    else if (z < 0) bulk1_count_dir--;
    double excess = std::max(0.0, std::abs(z) - tau);
    sum_excess += excess;
  }

  us.sum_bulk1 = total_b1;
  us.sum_bulk2 = total_b2;
  us.dominant_dir = (bulk1_count_dir >= 0) ? int8_t(+1) : int8_t(-1);
  int n_dominant = (us.dominant_dir > 0)
    ? (us.n_kmers + bulk1_count_dir) / 2
    : (us.n_kmers - bulk1_count_dir) / 2;
  us.coherence = static_cast<float>(n_dominant) / us.n_kmers;
  us.mean_z = static_cast<float>(sum_abs_z / us.n_kmers);
  us.sum_excess_mass = static_cast<float>(sum_excess);
  us.excess_per_kb = static_cast<float>(sum_excess / (us.length / 1000.0));
  us.rank_score = us.excess_per_kb * us.coherence;
  return us;
}

inline std::vector<UnitigScore> score_and_rank_unitigs_kmc(
    const std::vector<std::string>& unitigs, uint32_t k,
    KmcRandomAccess& db1, KmcRandomAccess& db2,
    const Params& params, float tau)
{
  std::vector<UnitigScore> scores;
  scores.reserve(unitigs.size());

  for (size_t i = 0; i < unitigs.size(); ++i) {
    auto us = score_unitig_kmc(unitigs[i], k, db1, db2, params, tau);
    us.id = "unitig_" + std::to_string(i + 1);
    us.sequence = unitigs[i];
    scores.push_back(std::move(us));
  }

  std::sort(scores.begin(), scores.end(), [](const UnitigScore& a, const UnitigScore& b) {
    if (a.rank_score != b.rank_score) return a.rank_score > b.rank_score;
    return a.sum_excess_mass > b.sum_excess_mass;
  });

  for (size_t i = 0; i < scores.size(); ++i) {
    const char* dir_str = (scores[i].dominant_dir > 0) ? "BULK1" : "BULK2";
    scores[i].id = "unitig_" + std::to_string(i + 1) + "_" + dir_str;
  }
  return scores;
}

}  // namespace kbsa

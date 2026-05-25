#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "../extern/kmc/kmc_api/kmc_file.h"

namespace kbsa {

struct Params
{
  double bulk1_to_bulk2_scale {1.0};
  uint64_t min_depth {5};
  uint64_t max_depth {500};
  double kai_min {0.8};
  double g_min {20.0};
  double error_rate {0.005};
};

struct HistEntry { uint64_t freq; uint64_t count; };

inline std::vector<HistEntry> load_histogram(const std::string& path)
{
  std::vector<HistEntry> hist;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    uint64_t f, c;
    if (iss >> f >> c) hist.push_back({f, c});
  }
  return hist;
}

inline std::vector<HistEntry> compute_histogram_from_db(const std::string& db_path)
{
  CKMCFile db;
  if (!db.OpenForListing(db_path))
    throw std::runtime_error("Cannot open KMC DB for histogram: " + db_path);

  CKMCFileInfo info;
  db.Info(info);
  CKmerAPI kmer(info.kmer_length);
  uint32 count = 0;

  std::vector<uint64_t> freq_counts;
  while (db.ReadNextKmer(kmer, count)) {
    if (count >= freq_counts.size())
      freq_counts.resize(count + 1, 0);
    freq_counts[count]++;
  }
  db.Close();

  std::vector<HistEntry> hist;
  for (size_t i = 1; i < freq_counts.size(); ++i) {
    if (freq_counts[i] > 0)
      hist.push_back({i, freq_counts[i]});
  }
  return hist;
}

inline uint64_t find_peak_depth(const std::vector<HistEntry>& hist)
{
  if (hist.empty()) return 0;

  // Find valley (local minimum) between error distribution and signal peak.
  // Error k-mers dominate low frequencies and decay; signal k-mers rise after.
  // Strategy: walk from freq=2 until count starts increasing → that's the valley.
  // Then find the maximum after the valley.

  // Sort by freq to ensure order
  std::vector<HistEntry> sorted_hist = hist;
  std::sort(sorted_hist.begin(), sorted_hist.end(),
    [](const HistEntry& a, const HistEntry& b) { return a.freq < b.freq; });

  // Find first entry with freq >= 2
  size_t start = 0;
  for (size_t i = 0; i < sorted_hist.size(); ++i) {
    if (sorted_hist[i].freq >= 2) { start = i; break; }
  }

  // Walk forward to find valley (first local minimum)
  uint64_t valley_freq = sorted_hist[start].freq;
  for (size_t i = start + 1; i < sorted_hist.size(); ++i) {
    if (sorted_hist[i].count > sorted_hist[i - 1].count) {
      valley_freq = sorted_hist[i - 1].freq;
      break;
    }
  }

  // Find peak after valley
  uint64_t peak_freq = valley_freq;
  uint64_t peak_count = 0;
  for (const auto& e : sorted_hist) {
    if (e.freq > valley_freq && e.count > peak_count) {
      peak_count = e.count;
      peak_freq = e.freq;
    }
  }

  // Fallback: if no valley found (monotonically increasing), use old behavior
  if (peak_count == 0) {
    for (const auto& e : hist) {
      if (e.freq >= 2 && e.count > peak_count) {
        peak_count = e.count;
        peak_freq = e.freq;
      }
    }
  }

  return peak_freq;
}

inline uint64_t total_kmers_above(const std::vector<HistEntry>& hist, uint64_t min_freq)
{
  uint64_t total = 0;
  for (const auto& e : hist)
    if (e.freq >= min_freq) total += e.count;
  return total;
}

inline Params auto_params(const std::string& bulk1_hist_path,
                          const std::string& bulk2_hist_path,
                          double kai_min = 0.8,
                          double g_min = 20.0,
                          double error_rate = 0.005)
{
  auto bulk1_hist = load_histogram(bulk1_hist_path);
  auto bulk2_hist = load_histogram(bulk2_hist_path);

  Params p;
  p.kai_min = kai_min;
  p.g_min = g_min;
  p.error_rate = error_rate;

  // Find peak depth (mode of true kmer distribution)
  uint64_t bulk1_peak = find_peak_depth(bulk1_hist);
  uint64_t bulk2_peak = find_peak_depth(bulk2_hist);
  uint64_t peak = std::max(bulk1_peak, bulk2_peak);

  // min_depth: adaptive cutoff (NIKS-style 2/3 of valley, minimum 3)
  // Simplified: use peak/5 as minimum, ensures low-count noise excluded
  p.min_depth = std::max<uint64_t>(3, peak / 5);

  // max_depth: 50x peak to exclude repetitive sequences
  p.max_depth = std::max<uint64_t>(p.min_depth * 2, peak * 50);

  // Scale factor: total_case / total_ctrl (depth normalization)
  uint64_t case_total = total_kmers_above(bulk1_hist, p.min_depth / 2);
  uint64_t ctrl_total = total_kmers_above(bulk2_hist, p.min_depth / 2);
  p.bulk1_to_bulk2_scale = ctrl_total > 0
    ? static_cast<double>(ctrl_total) / static_cast<double>(case_total)
    : 1.0;

  return p;
}

}  // namespace kbsa

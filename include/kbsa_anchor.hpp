#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include "kbsa_core.hpp"
#include "kbsa_params.hpp"
#include "../extern/kmc/kmc_api/kmc_file.h"

namespace kbsa {

struct AnchorHit
{
  uint64_t position;
  double kai_reg;
  double g_score;
  Significance direction;
};

struct CandidateInterval
{
  std::string chrom;
  uint64_t start;
  uint64_t end;
  uint32_t hit_count;
  uint32_t window_kmers;
  double density;
  double consistency;
  double mean_kai_deviation;
  double max_g_score;
};

inline void reverse_complement(const char* seq, char* out, uint32_t len)
{
  static const auto kComplement = []() {
    std::array<char, 256> t {};
    t['A'] = 'T'; t['T'] = 'A'; t['C'] = 'G'; t['G'] = 'C';
    t['a'] = 'T'; t['t'] = 'A'; t['c'] = 'G'; t['g'] = 'C';
    return t;
  }();
  for (uint32_t i = 0; i < len; ++i) {
    char c = kComplement[static_cast<unsigned char>(seq[len - 1 - i])];
    out[i] = c ? c : 'N';
  }
  out[len] = '\0';
}

inline void canonicalize(const char* fwd, char* canon, char* rc_buf, uint32_t k)
{
  reverse_complement(fwd, rc_buf, k);
  if (memcmp(fwd, rc_buf, k) <= 0)
    memcpy(canon, fwd, k + 1);
  else
    memcpy(canon, rc_buf, k + 1);
}

class KmcRandomAccess
{
public:
  KmcRandomAccess() = default;

  bool open(const std::string& path)
  {
    if (!m_db.OpenForRA(path)) return false;
    CKMCFileInfo info;
    m_db.Info(info);
    m_kmer_len = info.kmer_length;
    m_kmer = CKmerAPI(m_kmer_len);
    return true;
  }

  void close() { m_db.Close(); }

  uint32_t kmer_length() const { return m_kmer_len; }

  uint64_t query(const char* kmer_str)
  {
    m_kmer.from_string(kmer_str);
    uint32 count = 0;
    if (m_db.CheckKmer(m_kmer, count))
      return count;
    return 0;
  }

private:
  CKMCFile m_db;
  CKmerAPI m_kmer;
  uint32_t m_kmer_len {31};
};

struct AnchorParams
{
  uint32_t window_size {5000};
  uint32_t step_size {1000};
  double min_density {0.01};
  double min_consistency {0.7};
  double kai_min {0.55};
  double g_min {3.0};
  double error_rate {0.005};
};

class FastaReader
{
public:
  explicit FastaReader(const std::string& path)
    : m_in(path)
  {
    if (!m_in.good())
      throw std::runtime_error("Cannot open FASTA: " + path);
  }

  bool next_chromosome(std::string& name, std::string& seq)
  {
    name.clear();
    seq.clear();

    std::string line;
    if (!m_header_buffered) {
      while (std::getline(m_in, line)) {
        if (!line.empty() && line[0] == '>') {
          m_buffered_header = line;
          m_header_buffered = true;
          break;
        }
      }
    }

    if (!m_header_buffered) return false;

    auto space_pos = m_buffered_header.find_first_of(" \t", 1);
    name = m_buffered_header.substr(1, space_pos - 1);
    m_header_buffered = false;

    while (std::getline(m_in, line)) {
      if (!line.empty() && line[0] == '>') {
        m_buffered_header = line;
        m_header_buffered = true;
        break;
      }
      seq.append(line);
    }

    for (auto& c : seq) c = static_cast<char>(toupper(c));
    return true;
  }

private:
  std::ifstream m_in;
  std::string m_buffered_header;
  bool m_header_buffered {false};
};

inline bool has_ambiguous(const char* seq, uint32_t len)
{
  static const auto kIsACGT = []() {
    std::array<bool, 256> t {};
    t['A'] = true; t['C'] = true; t['G'] = true; t['T'] = true;
    return t;
  }();
  for (uint32_t i = 0; i < len; ++i) {
    if (!kIsACGT[static_cast<unsigned char>(seq[i])]) return true;
  }
  return false;
}

}  // namespace kbsa

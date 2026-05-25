#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>

#include "../extern/kmc/kmc_api/kmc_file.h"

namespace kbsa {

struct MergedKmer
{
  const char* kmer_str {nullptr};
  uint32_t kmer_len {0};
  uint64_t bulk1_count {0};
  uint64_t bulk2_count {0};
};

class KmcMergeIterator
{
public:
  KmcMergeIterator(const std::string& bulk1_db_path,
                   const std::string& bulk2_db_path)
  {
    if (!m_bulk1_db.OpenForListing(bulk1_db_path))
      throw std::runtime_error("Cannot open bulk1 db: " + bulk1_db_path);
    if (!m_bulk2_db.OpenForListing(bulk2_db_path))
      throw std::runtime_error("Cannot open bulk2 db: " + bulk2_db_path);

    CKMCFileInfo b1_info, b2_info;
    m_bulk1_db.Info(b1_info);
    m_bulk2_db.Info(b2_info);

    if (b1_info.kmer_length != b2_info.kmer_length)
      throw std::runtime_error("K-mer length mismatch: bulk1=" +
        std::to_string(b1_info.kmer_length) + " bulk2=" +
        std::to_string(b2_info.kmer_length));

    m_kmer_len = b1_info.kmer_length;
    m_bulk1_str.assign(m_kmer_len + 1, '\0');
    m_bulk2_str.assign(m_kmer_len + 1, '\0');
    m_bulk1_kmer = CKmerAPI(m_kmer_len);
    m_bulk2_kmer = CKmerAPI(m_kmer_len);

    m_bulk1_valid = m_bulk1_db.ReadNextKmer(m_bulk1_kmer, m_bulk1_count);
    m_bulk2_valid = m_bulk2_db.ReadNextKmer(m_bulk2_kmer, m_bulk2_count);
    if (m_bulk1_valid) m_bulk1_kmer.to_string(m_bulk1_str.data());
    if (m_bulk2_valid) m_bulk2_kmer.to_string(m_bulk2_str.data());
  }

  ~KmcMergeIterator()
  {
    m_bulk1_db.Close();
    m_bulk2_db.Close();
  }

  KmcMergeIterator(const KmcMergeIterator&) = delete;
  KmcMergeIterator& operator=(const KmcMergeIterator&) = delete;

  bool next(MergedKmer& out)
  {
    if (!m_bulk1_valid && !m_bulk2_valid) return false;

    if (!m_bulk2_valid) {
      out = {m_bulk1_str.data(), m_kmer_len, m_bulk1_count, 0};
      advance_bulk1();
      return true;
    }
    if (!m_bulk1_valid) {
      out = {m_bulk2_str.data(), m_kmer_len, 0, m_bulk2_count};
      advance_bulk2();
      return true;
    }

    int cmp = memcmp(m_bulk1_str.data(), m_bulk2_str.data(), m_kmer_len);
    if (cmp < 0) {
      out = {m_bulk1_str.data(), m_kmer_len, m_bulk1_count, 0};
      advance_bulk1();
    } else if (cmp > 0) {
      out = {m_bulk2_str.data(), m_kmer_len, 0, m_bulk2_count};
      advance_bulk2();
    } else {
      out = {m_bulk1_str.data(), m_kmer_len, m_bulk1_count, m_bulk2_count};
      advance_bulk1();
      advance_bulk2();
    }
    return true;
  }

  uint32_t kmer_length() const { return m_kmer_len; }

private:
  void advance_bulk1()
  {
    m_bulk1_valid = m_bulk1_db.ReadNextKmer(m_bulk1_kmer, m_bulk1_count);
    if (m_bulk1_valid) m_bulk1_kmer.to_string(m_bulk1_str.data());
  }

  void advance_bulk2()
  {
    m_bulk2_valid = m_bulk2_db.ReadNextKmer(m_bulk2_kmer, m_bulk2_count);
    if (m_bulk2_valid) m_bulk2_kmer.to_string(m_bulk2_str.data());
  }

  CKMCFile m_bulk1_db, m_bulk2_db;
  CKmerAPI m_bulk1_kmer, m_bulk2_kmer;
  uint32 m_bulk1_count {0}, m_bulk2_count {0};
  bool m_bulk1_valid {false}, m_bulk2_valid {false};
  uint32_t m_kmer_len {31};
  std::vector<char> m_bulk1_str, m_bulk2_str;
};

}  // namespace kbsa

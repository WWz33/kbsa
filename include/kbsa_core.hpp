#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace kbsa {

enum class Significance : uint8_t { NO = 0, BULK2 = 1, BULK1 = 2 };

inline const char* significance_str(Significance s) noexcept
{
  switch (s) {
    case Significance::BULK2: return "BULK2";
    case Significance::BULK1:    return "BULK1";
    default:                    return "NO";
  }
}

struct Score
{
  double kai_reg {0.5};
  double g_score {0.0};
  double z_score {0.0};
  Significance sign {Significance::NO};
  bool passed {false};
  uint64_t cnt_bulk2_raw {0};
  uint64_t cnt_bulk1_raw {0};
  double cnt_bulk2_adj {0.0};
  double cnt_bulk1_adj {0.0};
};

inline double dynamic_epsilon(uint64_t count, double error_rate) noexcept
{
  return std::max(1.0, static_cast<double>(count) * error_rate);
}

inline double regularized_kai(double bulk1_count,
                              double bulk2_count,
                              double bulk1_to_bulk2_scale,
                              double bulk1_epsilon,
                              double bulk2_epsilon) noexcept
{
  const double adjusted_bulk1 = bulk1_count * bulk1_to_bulk2_scale;
  const double denom = adjusted_bulk1 + bulk1_epsilon + bulk2_count + bulk2_epsilon;
  if (denom <= 0.0) return 0.5;
  return (adjusted_bulk1 + bulk1_epsilon) / denom;
}

inline double g_test(double adjusted_control, double adjusted_case) noexcept
{
  if (adjusted_control <= 0.0 || adjusted_case <= 0.0) return 0.0;
  const double total = adjusted_control + adjusted_case;
  if (total <= 0.0) return 0.0;
  const double expected = total / 2.0;
  if (expected <= 0.0) return 0.0;
  const double g = 2.0 * (
    adjusted_control * std::log(adjusted_control / expected) +
    adjusted_case * std::log(adjusted_case / expected)
  );
  if (!std::isfinite(g) || g < 0.0) return 0.0;
  return g;
}

inline Score score_kmer(uint64_t bulk2_count,
                        uint64_t bulk1_count,
                        double bulk1_to_bulk2_scale,
                        uint64_t min_depth,
                        uint64_t max_depth,
                        double kai_min,
                        double g_min,
                        double error_rate) noexcept
{
  Score s {};
  const auto total = bulk2_count + bulk1_count;
  if (total < min_depth || total > max_depth) return s;

  const double bulk1_eps = dynamic_epsilon(bulk1_count, error_rate);
  const double bulk2_eps = dynamic_epsilon(bulk2_count, error_rate);

  s.kai_reg = regularized_kai(
    static_cast<double>(bulk1_count),
    static_cast<double>(bulk2_count),
    bulk1_to_bulk2_scale, bulk1_eps, bulk2_eps);

  const double adj_bulk1 = static_cast<double>(bulk1_count) * bulk1_to_bulk2_scale + bulk1_eps;
  const double adj_bulk2 = static_cast<double>(bulk2_count) + bulk2_eps;
  s.g_score = g_test(adj_bulk1, adj_bulk2);
  s.cnt_bulk1_raw = bulk1_count;
  s.cnt_bulk2_raw = bulk2_count;
  s.cnt_bulk1_adj = adj_bulk1;
  s.cnt_bulk2_adj = adj_bulk2;

  // kai_reg = bulk1 / (bulk1 + bulk2): >0.5 = BULK1-enriched
  s.sign = (s.kai_reg >= 0.5) ? Significance::BULK1 : Significance::BULK2;

  // z > 0 = BULK1-enriched, z < 0 = BULK2-enriched
  const double sign_val = (s.sign == Significance::BULK1) ? 1.0 : -1.0;
  s.z_score = sign_val * std::sqrt(s.g_score);

  // Thresholds only gate emission (passed), never direction.
  const bool kai_pass = (s.kai_reg >= kai_min) || (s.kai_reg <= (1.0 - kai_min));
  s.passed = kai_pass && (s.g_score >= g_min);
  return s;
}

}  // namespace kbsa

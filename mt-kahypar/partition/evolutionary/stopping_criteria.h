#pragma once
#include "mt-kahypar/partition/evolutionary/evo_logs.h"
#include <vector>
#include <utility>

namespace mt_kahypar::evolutionary::stopping {

using mt_kahypar::evolutionary::IterationLogEntry;
using mt_kahypar::evolutionary::ImprovementLogEntry;

/**
 * Sliding-window improvement rate stopping criterion.
 *
 * Given a sequence of (iteration, objective) pairs (e.g., km1 values) and
 * parameters (early_window_improvs, recent_window_improvs, alpha,
 * max_iters_without_improv), this function decides when to stop.
 */

bool sliding_window_improvement_rate_stop(
  const int last_improvement_iteration,
  const std::vector<ImprovementLogEntry>& improvement_log_entries,
  const std::vector<IterationLogEntry>& iteration_log_entries,
  const int early_window_improvs,
  const int recent_window_improvs,
  const double alpha,
  const int max_iters_without_improv,
  double& early_window_improvement_rate);

}  // namespace mt_kahypar::evolutionary::stopping
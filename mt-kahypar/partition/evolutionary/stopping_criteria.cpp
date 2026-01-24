#include "mt-kahypar/partition/evolutionary/stopping_criteria.h"
#include <iostream>
#include <algorithm>
#include <limits>

namespace mt_kahypar::evolutionary::stopping {


bool sliding_window_improvement_rate_stop(
  const int last_improvement_iteration,
  const std::vector<ImprovementLogEntry>& improvement_log_entries,
  const std::vector<IterationLogEntry>& iteration_log_entries,
  const int early_window_improvs,
  const int recent_window_improvs,
  const double alpha,
  const int max_iters_without_improv,
  double& early_window_improvement_rate)
{
    if (iteration_log_entries.empty()) {
        return false;
    }

    const int last_iter = iteration_log_entries.back().iteration;
    
    // Plateau fallback: too many iterations since last improvement
    if (last_improvement_iteration >= 0 &&
        last_improvement_iteration + max_iters_without_improv < last_iter) {
        return true;
    }

    if (static_cast<int>(improvement_log_entries.size()) < early_window_improvs) {
        return false;
    }

    // Identify evolutionary (non-"Initial") improvements
    // Only needs to be done if early rate not yet computed
    std::vector<int> evo_indices;
    if (early_window_improvement_rate < 0.0) {
        for (int i = 0; i < static_cast<int>(improvement_log_entries.size()); ++i) {
            if (improvement_log_entries[i].operation_type != "Initial") {
                evo_indices.push_back(i);
            }
        }
    }

    if (evo_indices.size() < early_window_improvs) {
        return false;
    }

    // Compute early window improvement rate once and write it back
    if (early_window_improvement_rate < 0.0) {
        // First early_window_improvs improvements define the early window
        const auto& first = improvement_log_entries[evo_indices[0]];
        const auto& last  = improvement_log_entries[evo_indices[early_window_improvs - 1]];

        const double earliest_km1      = first.km1;
        const double latest_km1        = last.km1;
        const int earliest_iteration   = first.iteration;
        const int latest_iteration     = last.iteration;
        const int span                 = std::max(latest_iteration - earliest_iteration, 1);

        early_window_improvement_rate = (earliest_km1 - latest_km1) / static_cast<double>(span);

        // This should actually never trigger because at least 2 improvements (span >= 1) are needed
        if (early_window_improvement_rate <= 0.0) {
            return false;
        }
    }

    const double early_rate = early_window_improvement_rate;

    // Compute recent window improvement rate using only non-"Initial" entries
    if (static_cast<int>(evo_indices.size()) <= early_window_improvs)
        return false;

    const int n_evo = static_cast<int>(evo_indices.size());
    const int end_idx = evo_indices[n_evo - 1];
    const int start_idx = evo_indices[std::max(0, n_evo - recent_window_improvs)];

    const auto& win_start = improvement_log_entries[start_idx];
    const auto& win_end   = improvement_log_entries[end_idx];

    const double win_start_km1 = win_start.km1;
    const double win_end_km1 = win_end.km1;
    const int win_start_iter = win_start.iteration;
    const int win_end_iter = win_end.iteration;
    const int win_span = std::max(win_end_iter - win_start_iter, 1);
    const double recent_rate = (win_start_km1 - win_end_km1) / static_cast<double>(win_span);

    if (recent_rate < alpha * early_rate) {
        return true;
    }

    return false;

}

}  // namespace mt_kahypar::evolutionary::stopping
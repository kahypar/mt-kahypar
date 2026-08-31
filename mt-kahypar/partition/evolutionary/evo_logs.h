#pragma once

#include <string>

namespace mt_kahypar::evolutionary {

struct IterationLogEntry {
  int iteration;
  long long timestamp;
  double km1;
};

struct ImprovementLogEntry {
  long long timestamp;
  int iteration;
  double km1;
  std::string operation_type;
};

struct PopulationSummaryLogEntry {
  int iteration;
  long long timestamp;
  double best_km1;
  double median_km1;
  double mean_km1;
  double worst_km1;
  size_t min_distance;
  double median_distance;
  double mean_distance;
  size_t max_distance;
};

struct OffspringLogEntry {
  int iteration;
  long long timestamp;
  std::string operator_name;
  double parent_best_km1;
  double offspring_km1;
  double fitness_delta;
  size_t dist_to_parent;
  int64_t dist_to_replaced;
  size_t min_dist_pop;
  bool was_accepted;
  bool was_new_global_best;
  int replaced_slot;
};

}  // namespace mt_kahypar::evolutionary
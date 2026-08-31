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

}  // namespace mt_kahypar::evolutionary
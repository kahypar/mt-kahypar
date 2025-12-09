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

}  // namespace mt_kahypar::evolutionary
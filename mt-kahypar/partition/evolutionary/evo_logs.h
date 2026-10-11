#pragma once

#include <chrono>
#include <string>

namespace mt_kahypar::evolutionary {

// Milliseconds since epoch -- the timestamp format of all evolutionary log files
inline long long currentTimestampMs() {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

struct ImprovementLogEntry {
  long long timestamp;
  int iteration;
  double km1;
  std::string operation_type;
};

}  // namespace mt_kahypar::evolutionary

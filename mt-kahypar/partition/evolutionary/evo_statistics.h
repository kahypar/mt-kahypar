#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <vector>

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/evolutionary/evo_logs.h"

namespace mt_kahypar::evolutionary {

enum class EvoOperator : uint8_t {
  basic_combine,
  edge_frequency_combine,
  synthetic_parent_combine,
  vcycle_mutation,
  new_ip_mutation,
  COUNT
};
constexpr size_t kNumEvoOperators = static_cast<size_t>(EvoOperator::COUNT);

EvoOperator toEvoOperator(EvoCombineStrategy strategy);
EvoOperator toEvoOperator(EvoMutateStrategy strategy);
const char* operatorName(EvoOperator op);

struct OperatorStats {
  uint64_t attempts = 0;
  uint64_t accepted = 0;
};

// Everything that is known about one offspring after its insertion attempt
struct OffspringEvent {
  int iteration;
  long long timestamp_ms;
  EvoOperator op;
  HyperedgeWeight parent_best_fitness;
  HyperedgeWeight fitness;
  size_t dist_to_parent;
  bool was_accepted;
  int64_t dist_to_replaced;  // -1 if not accepted
  size_t min_dist_pop;
  int replaced_slot;  // -1 if not accepted
};

// Statistics of one evolutionary run.
class EvoStatistics {
 public:
  EvoStatistics();

  // Individuals of the initial population: only update the global best.
  // Returns true if the individual is a new global best.
  bool recordInitial(HyperedgeWeight fitness, int iteration, long long timestamp_ms);
  // Returns true if the offspring is a new global best.
  bool recordOffspring(const OffspringEvent& event);

  OperatorStats operatorStats(EvoOperator op) const; 
  // Every new global best (including the initial population) in chronological order
  std::vector<ImprovementLogEntry> improvements() const;

 private:
  mutable std::mutex _mutex;
  HyperedgeWeight _global_best;
  std::array<OperatorStats, kNumEvoOperators> _operators;  // indexed by EvoOperator
  std::vector<ImprovementLogEntry> _improvements;
};

}  // namespace mt_kahypar::evolutionary

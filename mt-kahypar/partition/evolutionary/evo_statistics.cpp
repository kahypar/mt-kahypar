#include "mt-kahypar/partition/evolutionary/evo_statistics.h"

#include <limits>

#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::evolutionary {

EvoOperator toEvoOperator(const EvoCombineStrategy strategy) {
  switch (strategy) {
    case EvoCombineStrategy::basic: return EvoOperator::basic_combine;
    case EvoCombineStrategy::edge_frequency: return EvoOperator::edge_frequency_combine;
    case EvoCombineStrategy::synthetic_parent: return EvoOperator::synthetic_parent_combine;
    default: throw UnsupportedOperationException("Undefined combine strategy");
  }
}

EvoOperator toEvoOperator(const EvoMutateStrategy strategy) {
  switch (strategy) {
    case EvoMutateStrategy::vcycle: return EvoOperator::vcycle_mutation;
    case EvoMutateStrategy::new_initial_partitioning_vcycle: return EvoOperator::new_ip_mutation;
    default: throw UnsupportedOperationException("Undefined mutation strategy");
  }
}

const char* operatorName(const EvoOperator op) {
  switch (op) {
    case EvoOperator::basic_combine: return "Basic_Combine";
    case EvoOperator::edge_frequency_combine: return "EdgeFreq_Combine";
    case EvoOperator::synthetic_parent_combine: return "Synthetic_Combine";
    case EvoOperator::vcycle_mutation: return "vCycle_Mutate";
    case EvoOperator::new_ip_mutation: return "NewInit_Mutate";
    default: return "Unknown";
  }
}

EvoStatistics::EvoStatistics() :
  _global_best(std::numeric_limits<HyperedgeWeight>::max()),
  _operators() { }

bool EvoStatistics::recordInitial(const HyperedgeWeight fitness, const int iteration, const long long timestamp_ms) {
  std::lock_guard<std::mutex> lock(_mutex);
  if (fitness < _global_best) {
    _global_best = fitness;
    _improvements.push_back({timestamp_ms, iteration, static_cast<double>(fitness), "Initial"});
    return true;
  }
  return false;
}

bool EvoStatistics::recordOffspring(const OffspringEvent& event) {
  std::lock_guard<std::mutex> lock(_mutex);
  OperatorStats& stats = _operators[static_cast<size_t>(event.op)];
  ++stats.attempts;
  if (event.was_accepted) ++stats.accepted;
  if (event.fitness < _global_best) {
    _global_best = event.fitness;
    _improvements.push_back({event.timestamp_ms, event.iteration, static_cast<double>(event.fitness), operatorName(event.op)});
    return true;
  }
  return false;
}

OperatorStats EvoStatistics::operatorStats(const EvoOperator op) const {
  std::lock_guard<std::mutex> lock(_mutex);
  return _operators[static_cast<size_t>(op)];
}

std::vector<ImprovementLogEntry> EvoStatistics::improvements() const {
  std::lock_guard<std::mutex> lock(_mutex);
  return _improvements;
}

}  // namespace mt_kahypar::evolutionary

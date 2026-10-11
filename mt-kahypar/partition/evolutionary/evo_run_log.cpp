#include "mt-kahypar/partition/evolutionary/evo_run_log.h"

#include <algorithm>
#include <fstream>
#include <iomanip>

namespace mt_kahypar::evolutionary {

EvoRunLog::EvoRunLog(const EvolutionaryParameters& params, const long long start_timestamp_ms) :
  _level(params.log_level),
  _interval(std::max<size_t>(params.log_interval, 1)),
  _start_timestamp_ms(start_timestamp_ms),
  _history_file(params.history_file),
  _offspring_log_file(params.offspring_log_file),
  _diff_matrix_file(params.diff_matrix_file),
  _population_summary_file(params.population_summary_file) { }

bool EvoRunLog::logsOffspring() const {
  return _level != EvoLogLevel::none && !_offspring_log_file.empty();
}

void EvoRunLog::recordOffspring(const OffspringEvent& event, const bool was_new_best) {
  std::lock_guard<std::mutex> lock(_mutex);
  _offspring.push_back({event, was_new_best});
}

bool EvoRunLog::logsPopulation() const {
  return _level == EvoLogLevel::full && (!_diff_matrix_file.empty() || !_population_summary_file.empty());
}

bool EvoRunLog::isPopulationSnapshotDue(const int iteration) const {
  return logsPopulation() && iteration % _interval == 0;
}

void EvoRunLog::recordDiffMatrix(const Population& population) {
  if (_level != EvoLogLevel::full || _diff_matrix_file.empty()) return;
  DiffMatrix diff_matrix = population.updateDiffMatrix();
  std::lock_guard<std::mutex> lock(_mutex);
  _diff_matrices.push_back(std::move(diff_matrix));
}

void EvoRunLog::recordPopulationSummary(const Population& population, const int iteration, const long long timestamp_ms) {
  if (_level != EvoLogLevel::full || _population_summary_file.empty()) return;
  const PopulationSummary summary = population.computeSummary();
  std::lock_guard<std::mutex> lock(_mutex);
  _population_snapshots.push_back({iteration, timestamp_ms, summary});
}

void EvoRunLog::writeFiles(const EvoStatistics& stats) const {
  if (_level == EvoLogLevel::none) return;
  std::lock_guard<std::mutex> lock(_mutex);
  writeHistory(stats);
  writeOffspringLog();
  writeDiffMatrices();
  writePopulationSummary();
}

// Same format as the history file of KaHyPar-E
void EvoRunLog::writeHistory(const EvoStatistics& stats) const {
  if (_history_file.empty()) return;
  std::ofstream out(_history_file, std::ios::out | std::ios::trunc);
  out << "Starttime: " << _start_timestamp_ms << "\n";
  for (const ImprovementLogEntry& entry : stats.improvements()) {
    out << entry.timestamp << ", " << entry.operation_type << ", " << static_cast<int64_t>(entry.km1) << "\n";
  }
}

void EvoRunLog::writeOffspringLog() const {
  if (!logsOffspring()) return;
  std::ofstream out(_offspring_log_file, std::ios::out | std::ios::trunc);
  out << "iteration,timestamp_ms,operator,parent_best_km1,offspring_km1,fitness_delta,dist_to_parent,"
         "dist_to_replaced,min_dist_pop,was_accepted,was_new_best,replaced_slot\n";
  for (const OffspringLogEntry& entry : _offspring) {
    const OffspringEvent& e = entry.event;
    out << e.iteration << ","
        << e.timestamp_ms << ","
        << operatorName(e.op) << ","
        << e.parent_best_fitness << ","
        << e.fitness << ","
        << (static_cast<int64_t>(e.parent_best_fitness) - static_cast<int64_t>(e.fitness)) << ","
        << e.dist_to_parent << ","
        << e.dist_to_replaced << ","
        << e.min_dist_pop << ","
        << (e.was_accepted ? 1 : 0) << ","
        << (entry.was_new_best ? 1 : 0) << ","
        << e.replaced_slot << "\n";
  }
}

// Same format as the diff matrix file of KaHyPar-E
void EvoRunLog::writeDiffMatrices() const {
  if (_level != EvoLogLevel::full || _diff_matrix_file.empty()) return;
  std::ofstream out(_diff_matrix_file, std::ios::out | std::ios::trunc);
  for (const DiffMatrix& diff_matrix : _diff_matrices) {
    for (const std::vector<size_t>& row : diff_matrix) {
      out << Population::toString(row) << "\n";
    }
    out << "---\n";
  }
}

void EvoRunLog::writePopulationSummary() const {
  if (_level != EvoLogLevel::full || _population_summary_file.empty()) return;
  std::ofstream out(_population_summary_file, std::ios::out | std::ios::trunc);
  out << "iteration,timestamp_ms,best_km1,median_km1,mean_km1,worst_km1,min_dist,median_dist,mean_dist,max_dist\n";
  out << std::fixed << std::setprecision(2);
  for (const PopulationSnapshot& snapshot : _population_snapshots) {
    const PopulationSummary& s = snapshot.summary;
    out << snapshot.iteration << "," << snapshot.timestamp_ms << ","
        << s.best_fitness << "," << s.median_fitness << ","
        << s.mean_fitness << "," << s.worst_fitness << ","
        << s.min_distance << "," << s.median_distance << ","
        << s.mean_distance << "," << s.max_distance << "\n";
  }
}

}  // namespace mt_kahypar::evolutionary

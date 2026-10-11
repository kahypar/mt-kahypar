#pragma once

#include <mutex>
#include <string>
#include <vector>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/evolutionary/evo_statistics.h"
#include "mt-kahypar/partition/evolutionary/population.h"

namespace mt_kahypar::evolutionary {

// Telemetry of one evolutionary run (for later analysis, never read by the algorithm).
// What is recorded depends on the log level and on which output files are set.
// Thread-safe; one instance per run (not static).
class EvoRunLog {
  struct OffspringLogEntry {
    OffspringEvent event;
    bool was_new_best;
  };

  struct PopulationSnapshot {
    int iteration;
    long long timestamp_ms;
    PopulationSummary summary;
  };

 public:
  EvoRunLog(const EvolutionaryParameters& params, long long start_timestamp_ms);

  // Offspring log (log level light or full)
  bool logsOffspring() const;
  void recordOffspring(const OffspringEvent& event, bool was_new_best);

  // Diff matrix and population summary (log level full)
  bool logsPopulation() const;
  // True if a population snapshot should be taken after the insertion of this iteration
  bool isPopulationSnapshotDue(int iteration) const;
  void recordDiffMatrix(const Population& population);
  void recordPopulationSummary(const Population& population, int iteration, long long timestamp_ms);

  // Writes all enabled files. The history file is built from the improvements in `stats`.
  void writeFiles(const EvoStatistics& stats) const;

 private:
  void writeHistory(const EvoStatistics& stats) const;
  void writeOffspringLog() const;
  void writeDiffMatrices() const;
  void writePopulationSummary() const;

  const EvoLogLevel _level;
  const size_t _interval;
  const long long _start_timestamp_ms;
  const std::string _history_file;
  const std::string _offspring_log_file;
  const std::string _diff_matrix_file;
  const std::string _population_summary_file;

  mutable std::mutex _mutex;
  std::vector<OffspringLogEntry> _offspring;
  std::vector<DiffMatrix> _diff_matrices;
  std::vector<PopulationSnapshot> _population_snapshots;
};

}  // namespace mt_kahypar::evolutionary

/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <mutex>
#include <string>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"

namespace mt_kahypar {
namespace utils {

struct InitialPartitionerSummary {

  explicit InitialPartitionerSummary(const InitialPartitioningAlgorithm algo) :
    algorithm(algo),
    total_sum_quality(0),
    total_time(0.0),
    total_best(0),
    total_calls(0) { }

  friend std::ostream & operator<< (std::ostream& str, const InitialPartitionerSummary& summary);

  void add(const InitialPartitionerSummary& summary) {
    ASSERT(algorithm == summary.algorithm);
    total_sum_quality += summary.total_sum_quality;
    total_time += summary.total_time;
    total_calls += summary.total_calls;
  }

  double average_quality() const {
    return static_cast<double>(total_sum_quality) / std::max(total_calls, 1UL);
  }

  double average_running_time() const {
    return static_cast<double>(total_time) / std::max(total_calls, 1UL);
  }

  double percentage_best(const size_t total_ip_calls) const {
    return ( static_cast<double>(total_best) / total_ip_calls ) * 100.0;
  }

  InitialPartitioningAlgorithm algorithm;
  double total_sum_quality;
  double total_time;
  size_t total_best;
  size_t total_calls;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitionerSummary& summary) {
  str << " avg_quality_" << summary.algorithm << "=" << summary.average_quality()
      << " total_time_" << summary.algorithm << "=" << summary.total_time
      << " total_best_" << summary.algorithm << "=" << summary.total_best;
  return str;
}

class InitialPartitioningStats {

 public:
  InitialPartitioningStats(const InitialPartitioningStats&) = delete;
  InitialPartitioningStats & operator= (const InitialPartitioningStats &) = delete;

  InitialPartitioningStats(InitialPartitioningStats&&) = delete;
  InitialPartitioningStats & operator= (InitialPartitioningStats &&) = delete;

  static InitialPartitioningStats& instance() {
    static InitialPartitioningStats instance;
    return instance;
  }

  void add_initial_partitioning_result(const InitialPartitioningAlgorithm best_algorithm,
                                       const size_t number_of_threads,
                                       const parallel::scalable_vector<InitialPartitionerSummary>& summary) {
    std::lock_guard<std::mutex> lock(_stat_mutex);
    ASSERT(summary.size() == _ip_summary.size());
    uint8_t best_algorithm_index = static_cast<uint8_t>(best_algorithm);
    ++_ip_summary[best_algorithm_index].total_best;

    for ( size_t i = 0; i < _num_initial_partitioner; ++i ) {
      _ip_summary[i].add(summary[i]);
    }

    _total_sum_number_of_threads += number_of_threads;
    ++_total_ip_calls;
  }

  double average_number_of_threads_per_ip_call() const {
    return static_cast<double>(_total_sum_number_of_threads) / _total_ip_calls;
  }

  void printInitialPartitioningStats() {
    LOG << "Initial Partitioning Algorithm Summary:";
    LOG << "Number of Initial Partitioning Calls =" << _total_ip_calls;
    LOG << "Average Number of Thread per IP Call ="
        << average_number_of_threads_per_ip_call() << "\n";
    std::cout << "\033[1m"
              << std::left << std::setw(30) << "Algorithm"
              << std::left << std::setw(15) << " Avg. Quality"
              << std::left << std::setw(15) << "  Total Time (s)"
              << std::left << std::setw(10) << "  Total Best"
              << std::left << std::setw(15) << " Total Best (%)"
              << "\033[0m" << std::endl;
    for ( const InitialPartitionerSummary& summary : _ip_summary ) {
      LOG << std::left << std::setw(30) << summary.algorithm
          << std::left << std::setw(15) << summary.average_quality()
          << std::left << std::setw(15) << summary.total_time
          << std::left << std::setw(10) << summary.total_best
          << std::left << std::setw(15) << summary.percentage_best(_total_ip_calls);
    }
  }

  friend std::ostream & operator<< (std::ostream& str, const InitialPartitioningStats& stats);

 private:
  explicit InitialPartitioningStats() :
    _stat_mutex(),
    _num_initial_partitioner(static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED)),
    _ip_summary(),
    _total_ip_calls(0),
    _total_sum_number_of_threads(0) {
    for ( uint8_t algo = 0; algo < _num_initial_partitioner; ++algo ) {
      _ip_summary.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
    }
  }

  std::mutex _stat_mutex;
  const uint8_t _num_initial_partitioner;
  parallel::scalable_vector<InitialPartitionerSummary> _ip_summary;
  size_t _total_ip_calls;
  size_t _total_sum_number_of_threads;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitioningStats& stats) {
  str << " average_number_of_threads_per_ip_call="
      << stats.average_number_of_threads_per_ip_call();
  for ( const InitialPartitionerSummary& summary : stats._ip_summary ) {
    str << summary;
  }
  return str;
}
}  // namespace utils
}  // namespace mt_kahypar

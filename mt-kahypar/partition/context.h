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

#include "kahypar/partition/context_enum_classes.h"
#include "kahypar/definitions.h"

namespace kahypar {

using namespace kahypar;

struct PartitioningParameters {
  Mode mode = Mode::UNDEFINED;
  Objective objective = Objective::UNDEFINED;
  double epsilon = std::numeric_limits<double>::max();
  PartitionID k = std::numeric_limits<PartitionID>::max();
  int seed = 0;  

  int time_limit = 0;
  HyperedgeID hyperedge_size_threshold = std::numeric_limits<HypernodeID>::max();

  bool verbose_output = false;
  bool quiet_mode = false;
  bool sp_process_output = false;
  bool write_partition_file = true;

  std::string graph_filename { };
  std::string graph_partition_filename { };
};

inline std::ostream& operator<< (std::ostream& str, const PartitioningParameters& params) {
  str << "Partitioning Parameters:" << std::endl;
  str << "  Hypergraph:                         " << params.graph_filename << std::endl;
  str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
  str << "  Mode:                               " << params.mode << std::endl;
  str << "  Objective:                          " << params.objective << std::endl;
  str << "  k:                                  " << params.k << std::endl;
  str << "  epsilon:                            " << params.epsilon << std::endl;
  str << "  seed:                               " << params.seed << std::endl;
  str << "  time limit:                         " << params.time_limit << "s" << std::endl;
  str << "  hyperedge size threshold:           " << params.hyperedge_size_threshold << std::endl;
  return str;
}

struct SharedMemoryParameters {
  size_t num_threads = 1;
};

inline std::ostream& operator<< (std::ostream& str, const SharedMemoryParameters& params) {
  str << "Shared Memory Parameters:             " << std::endl;
  str << "  Number of Threads:                  " << params.num_threads << std::endl;
  return str;
}

class Context {
 public:
  PartitioningParameters partition { };
  SharedMemoryParameters shared_memory { };
  ContextType type = ContextType::main;

  Context() { }

  ~Context() = default;

  Context(const Context& other) :
    partition(other.partition),
    type(other.type) { }

  Context& operator= (const Context&) = delete;

  bool isMainRecursiveBisection() const {
    return partition.mode == Mode::recursive_bisection && type == ContextType::main;
  }
};

inline std::ostream& operator<< (std::ostream& str, const Context& context) {
  str << "*******************************************************************************\n"
      << "*                            Partitioning Context                             *\n"
      << "*******************************************************************************\n"
      << context.partition
      << "-------------------------------------------------------------------------------\n"
      << context.shared_memory
      << "-------------------------------------------------------------------------------";
  return str;
}

} // namespace kahypar
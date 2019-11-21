/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#include "mt-kahypar/partition/context_enum_classes.h"
#include "kahypar/definitions.h"

namespace mt_kahypar {

using namespace kahypar;

struct PartitioningParameters {
  Mode mode = Mode::UNDEFINED;
  Objective objective = Objective::UNDEFINED;
  double epsilon = std::numeric_limits<double>::max();
  PartitionID k = std::numeric_limits<PartitionID>::max();
  int seed = 0;

  int time_limit = 0;
  std::vector<HypernodeWeight> perfect_balance_part_weights;
  std::vector<HypernodeWeight> max_part_weights;
  HyperedgeID hyperedge_size_threshold = 1000;

  bool verbose_output = false;
  bool quiet_mode = false;
  bool detailed_timings = false;
  bool sp_process_output = false;
  bool write_partition_file = false;

  std::string graph_filename { };
  std::string graph_partition_filename { };
  std::string graph_community_filename { };
};

inline std::ostream& operator<< (std::ostream& str, const PartitioningParameters& params) {
  str << "Partitioning Parameters:" << std::endl;
  str << "  Hypergraph:                         " << params.graph_filename << std::endl;
  str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
  str << "  Community File:                     " << params.graph_community_filename << std::endl;
  str << "  Mode:                               " << params.mode << std::endl;
  str << "  Objective:                          " << params.objective << std::endl;
  str << "  k:                                  " << params.k << std::endl;
  str << "  epsilon:                            " << params.epsilon << std::endl;
  str << "  seed:                               " << params.seed << std::endl;
  str << "  time limit:                         " << params.time_limit << "s" << std::endl;
  str << "  hyperedge size threshold:           " << params.hyperedge_size_threshold << std::endl;
  return str;
}

struct RatingParameters {
  RatingFunction rating_function = RatingFunction::UNDEFINED;
  HeavyNodePenaltyPolicy heavy_node_penalty_policy = HeavyNodePenaltyPolicy::UNDEFINED;
  AcceptancePolicy acceptance_policy = AcceptancePolicy::UNDEFINED;
};

inline std::ostream& operator<< (std::ostream& str, const RatingParameters& params) {
  str << "  Rating Parameters:" << std::endl;
  str << "    Rating Function:                  " << params.rating_function << std::endl;
  str << "    Heavy Node Penalty:               " << params.heavy_node_penalty_policy << std::endl;
  str << "    Acceptance Policy:                " << params.acceptance_policy << std::endl;
  return str;
}

struct PreprocessingParameters {
	CommunityDetectionStarExpansionWeightModification edge_weight_modification = CommunityDetectionStarExpansionWeightModification::hybrid;
};

struct CoarseningParameters {
  CoarseningAlgorithm algorithm = CoarseningAlgorithm::UNDEFINED;
  RatingParameters rating = { };
  HypernodeID contraction_limit_multiplier = std::numeric_limits<HypernodeID>::max();
  double max_allowed_weight_multiplier = std::numeric_limits<double>::max();
  bool use_hypernode_degree_threshold = false;

  // Those will be determined dynamically
  HypernodeWeight max_allowed_node_weight = 0;
  HypernodeID contraction_limit = 0;
  double hypernode_weight_fraction = 0.0;
  HyperedgeID hypernode_degree_threshold = std::numeric_limits<HyperedgeID>::max();
};


inline std::ostream& operator<< (std::ostream& str, const CoarseningParameters& params) {
  str << "Coarsening Parameters:" << std::endl;
  str << "  Algorithm:                          " << params.algorithm << std::endl;
  str << "  max allowed weight multiplier:      " << params.max_allowed_weight_multiplier << std::endl;
  str << "  contraction limit multiplier:       " << params.contraction_limit_multiplier << std::endl;
  str << "  use hypernode degree threshold:     " << std::boolalpha << params.use_hypernode_degree_threshold << std::endl;
  str << std::endl << params.rating;
  return str;
}

struct InitialPartitioningParameters {
  std::string context_file = "";
  bool call_kahypar_multiple_times = false;
  size_t runs = 1;
};

inline std::ostream& operator<< (std::ostream& str, const InitialPartitioningParameters& params) {
  str << "Initial Partitioning Parameters:" << std::endl;
  str << "  Initial Partitioning Context:       " << params.context_file << std::endl;
  str << "  Call KaHyPar multiple times:        " << std::boolalpha << params.call_kahypar_multiple_times << std::endl;
  str << "  Number of Runs:                     " << params.runs << std::endl;
  return str;
}

struct LabelPropagationParameters {
  LabelPropagationAlgorithm algorithm = LabelPropagationAlgorithm::do_nothing;
  size_t maximum_iterations = 1;
  size_t part_weight_update_frequency = 100;
  bool numa_aware = false;
  bool rebalancing = true;
  ExecutionType execution_policy = ExecutionType::UNDEFINED;
  double execution_policy_alpha = 2.0;
};

inline std::ostream& operator<< (std::ostream& str, const LabelPropagationParameters& params) {
  str << "  Label Propagation Parameters:" << std::endl;
  str << "    Algorithm:                        " << params.algorithm << std::endl;
  str << "    Maximum Iterations:               " << params.maximum_iterations << std::endl;
  str << "    Part Weight Update Frequency:     " << params.part_weight_update_frequency << std::endl;
  str << "    Numa Aware:                       " << std::boolalpha << params.numa_aware << std::endl;
  str << "    Rebalancing:                      " << std::boolalpha << params.rebalancing << std::endl;
  str << "    Execution Policy:                 " << params.execution_policy << std::endl;
  str << "    Execution Policy Alpha:           " << params.execution_policy_alpha << std::endl;
  return str;
}

struct RefinementParameters {
  LabelPropagationParameters label_propagation;
  bool use_batch_uncontractions = false;
  size_t batch_size = 1000;
};

inline std::ostream& operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  str << "  Use Batch Uncontractions:           " << std::boolalpha << params.use_batch_uncontractions << std::endl;
  if ( params.use_batch_uncontractions ) {
    str << "  Batch Size:                         " << params.batch_size << std::endl;
  }
  str << std::endl << params.label_propagation;
  return str;
}

struct SharedMemoryParameters {
  size_t num_threads = 1;
  bool use_community_redistribution = false;
  InitialHyperedgeDistribution initial_distribution = InitialHyperedgeDistribution::UNDEFINED;
  CommunityAssignmentObjective assignment_objective = CommunityAssignmentObjective::UNDEFINED;
  CommunityAssignmentStrategy assignment_strategy = CommunityAssignmentStrategy::UNDEFINED;
};

inline std::ostream& operator<< (std::ostream& str, const SharedMemoryParameters& params) {
  str << "Shared Memory Parameters:             " << std::endl;
  str << "  Number of Threads:                  " << params.num_threads << std::endl;
  str << "  Use Community Redistribution:       " << std::boolalpha << params.use_community_redistribution << std::endl;
  str << "  Initial Hyperedge Distribution:     " << params.initial_distribution << std::endl;
  str << "  Community Assignment Objective:     " << params.assignment_objective << std::endl;
  str << "  Community Assignment Strategy:      " << params.assignment_strategy << std::endl;
  return str;
}

class Context {
 public:
  PartitioningParameters partition { };
  PreprocessingParameters preprocessing { };
  CoarseningParameters coarsening { };
  InitialPartitioningParameters initial_partitioning { };
  RefinementParameters refinement { };
  SharedMemoryParameters shared_memory { };
  ContextType type = ContextType::main;

  Context() { }

  Context(const Context& other) :
    partition(other.partition),
    type(other.type) { }

  Context& operator= (const Context&) = delete;

  bool isMainRecursiveBisection() const {
    return partition.mode == Mode::recursive_bisection && type == ContextType::main;
  }

  void setupPartWeights(const HypernodeWeight total_hypergraph_weight) {
    partition.perfect_balance_part_weights.clear();
    partition.perfect_balance_part_weights.push_back(ceil(
                                                        total_hypergraph_weight
                                                        / static_cast<double>(partition.k)));
    for (PartitionID part = 1; part != partition.k; ++part) {
      partition.perfect_balance_part_weights.push_back(
        partition.perfect_balance_part_weights[0]);
    }
    partition.max_part_weights.clear();
    partition.max_part_weights.push_back((1 + partition.epsilon)
                                          * partition.perfect_balance_part_weights[0]);
    for (PartitionID part = 1; part != partition.k; ++part) {
      partition.max_part_weights.push_back(partition.max_part_weights[0]);
    }
  }
};

inline std::ostream& operator<< (std::ostream& str, const Context& context) {
  str << "*******************************************************************************\n"
      << "*                            Partitioning Context                             *\n"
      << "*******************************************************************************\n"
      << context.partition
      << "-------------------------------------------------------------------------------\n"
      << context.coarsening
      << "-------------------------------------------------------------------------------\n"
      << context.initial_partitioning
      << "-------------------------------------------------------------------------------\n"
      << context.refinement
      << "-------------------------------------------------------------------------------\n"
      << context.shared_memory
      << "-------------------------------------------------------------------------------";
  return str;
}

} // namespace mt_kahypar
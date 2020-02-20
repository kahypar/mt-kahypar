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

#include "kahypar/definitions.h"
#include "kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"

namespace mt_kahypar {
struct PartitioningParameters {
  Paradigm paradigm = Paradigm::nlevel;
  kahypar::Mode mode = kahypar::Mode::UNDEFINED;
  kahypar::Objective objective = kahypar::Objective::UNDEFINED;
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
  bool enable_progress_bar = false;
  bool sp_process_output = false;
  bool write_partition_file = false;

  std::string graph_filename { };
  std::string graph_partition_filename { };
  std::string graph_community_filename { };
};

inline std::ostream & operator<< (std::ostream& str, const PartitioningParameters& params) {
  str << "Partitioning Parameters:" << std::endl;
  str << "  Hypergraph:                         " << params.graph_filename << std::endl;
  str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
  str << "  Community File:                     " << params.graph_community_filename << std::endl;
  str << "  Paradigm:                           " << params.paradigm << std::endl;
  str << "  Mode:                               " << params.mode << std::endl;
  str << "  Objective:                          " << params.objective << std::endl;
  str << "  k:                                  " << params.k << std::endl;
  str << "  epsilon:                            " << params.epsilon << std::endl;
  str << "  seed:                               " << params.seed << std::endl;
  str << "  time limit:                         " << params.time_limit << "s" << std::endl;
  str << "  hyperedge size threshold:           " << params.hyperedge_size_threshold << std::endl;
  return str;
}

struct CommunityDetectionParameters {
  CommunityLoadBalancingStrategy load_balancing_strategy = CommunityLoadBalancingStrategy::none;
  size_t size_constraint_factor = 0;
  LouvainEdgeWeight edge_weight_function = LouvainEdgeWeight::UNDEFINED;
  uint32_t max_pass_iterations = std::numeric_limits<uint32_t>::max();
  long double min_eps_improvement = std::numeric_limits<long double>::max();
};

inline std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Load Balancing Strategy:          " << params.load_balancing_strategy << std::endl;
  if (params.load_balancing_strategy == CommunityLoadBalancingStrategy::size_constraint) {
    str << "    Size Constraint Factor:           " << params.size_constraint_factor << std::endl;
  }
  str << "    Edge Weight Function:             " << params.edge_weight_function << std::endl;
  str << "    Maximum Louvain-Pass Iterations:  " << params.max_pass_iterations << std::endl;
  str << "    Minimum Quality Improvement:      " << params.min_eps_improvement << std::endl;
  return str;
}

struct CommunityRedistributionParameters {
  CommunityAssignmentObjective assignment_objective = CommunityAssignmentObjective::UNDEFINED;
  CommunityAssignmentStrategy assignment_strategy = CommunityAssignmentStrategy::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const CommunityRedistributionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Community Assignment Objective:   " << params.assignment_objective << std::endl;
  str << "    Community Assignment Strategy:    " << params.assignment_strategy << std::endl;
  return str;
}

struct PreprocessingParameters {
  bool use_community_detection = false;
  bool use_community_redistribution = false;
  bool use_community_structure_from_file = false;
  CommunityDetectionParameters community_detection = { };
  CommunityRedistributionParameters community_redistribution = { };
};

inline std::ostream & operator<< (std::ostream& str, const PreprocessingParameters& params) {
  str << "Preprocessing Parameters:" << std::endl;
  str << "  Use Community Detection:            " << std::boolalpha << params.use_community_detection << std::endl;
  str << "  Use Community Redistribution:       " << std::boolalpha << params.use_community_redistribution << std::endl;
  str << "  Use Community Structure from File:  " << std::boolalpha << params.use_community_structure_from_file << std::endl;
  if (!params.use_community_structure_from_file && params.use_community_detection) {
    str << std::endl << params.community_detection;
  }
  if ( params.use_community_redistribution ) {
    str << std::endl << params.community_redistribution;
  }
  return str;
}

struct RatingParameters {
  RatingFunction rating_function = RatingFunction::UNDEFINED;
  HeavyNodePenaltyPolicy heavy_node_penalty_policy = HeavyNodePenaltyPolicy::UNDEFINED;
  AcceptancePolicy acceptance_policy = AcceptancePolicy::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const RatingParameters& params) {
  str << "  Rating Parameters:" << std::endl;
  str << "    Rating Function:                  " << params.rating_function << std::endl;
  str << "    Heavy Node Penalty:               " << params.heavy_node_penalty_policy << std::endl;
  str << "    Acceptance Policy:                " << params.acceptance_policy << std::endl;
  return str;
}

struct CoarseningParameters {
  CoarseningAlgorithm algorithm = CoarseningAlgorithm::UNDEFINED;
  RatingParameters rating = { };
  HypernodeID contraction_limit_multiplier = std::numeric_limits<HypernodeID>::max();
  double max_allowed_weight_multiplier = std::numeric_limits<double>::max();
  double max_allowed_high_degree_node_weight_multiplier = std::numeric_limits<double>::max();
  double multilevel_shrink_factor = std::numeric_limits<double>::max();
  bool ignore_already_matched_vertices = false;
  bool use_high_degree_vertex_threshold = false;

  // Those will be determined dynamically
  HypernodeWeight max_allowed_node_weight = 0;
  HypernodeWeight max_allowed_high_degree_node_weight = 0;
  HypernodeID contraction_limit = 0;
  HyperedgeID high_degree_vertex_threshold = std::numeric_limits<HyperedgeID>::max();
};

inline std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params) {
  str << "Coarsening Parameters:" << std::endl;
  str << "  Algorithm:                          " << params.algorithm << std::endl;
  str << "  max allowed weight multiplier:      " << params.max_allowed_weight_multiplier << std::endl;
  str << "  max allowed high degree multiplier: " << params.max_allowed_high_degree_node_weight_multiplier << std::endl;
  str << "  maximum allowed hypernode weight:   " << params.max_allowed_node_weight << std::endl;
  str << "  maximum allowed high-degree weight: " << params.max_allowed_high_degree_node_weight << std::endl;
  str << "  contraction limit multiplier:       " << params.contraction_limit_multiplier << std::endl;
  str << "  contraction limit:                  " << params.contraction_limit << std::endl;
  if ( params.algorithm == CoarseningAlgorithm::multilevel_coarsener ) {
    str << "  multilevel shrink factor:           " << params.multilevel_shrink_factor << std::endl;
    str << "  ignore already matched vertices:    " << std::boolalpha << params.ignore_already_matched_vertices << std::endl;
  }
  if ( params.use_high_degree_vertex_threshold ) {
    str << "  high degree vertex threshold:       " << params.high_degree_vertex_threshold << std::endl;
  }
  str << std::endl << params.rating;
  return str;
}

struct InitialPartitioningParameters {
  InitialPartitioningMode mode = InitialPartitioningMode::UNDEFINED;
  size_t runs = 1;
  bool use_adaptive_epsilon = false;
  size_t lp_maximum_iterations = 1;
  size_t lp_initial_block_size = 1;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params) {
  str << "Initial Partitioning Parameters:" << std::endl;
  str << "  Initial Partitioning Mode:          " << params.mode << std::endl;
  str << "  Number of Runs:                     " << params.runs << std::endl;
  str << "  Use Adaptive Epsilon:               " << std::boolalpha << params.use_adaptive_epsilon << std::endl;
  str << "  Maximum Iterations of LP IP:        " << params.lp_maximum_iterations << std::endl;
  str << "  Initial Block Size of LP IP:        " << params.lp_initial_block_size << std::endl;
  return str;
}

struct LabelPropagationParameters {
  LabelPropagationAlgorithm algorithm = LabelPropagationAlgorithm::do_nothing;
  size_t maximum_iterations = 1;
  double part_weight_update_factor = 0.01;
  bool localized = false;
  bool numa_aware = false;
  bool rebalancing = true;
  ExecutionType execution_policy = ExecutionType::UNDEFINED;
  double execution_policy_alpha = 2.0;
  bool execute_sequential = false;
};

inline std::ostream & operator<< (std::ostream& str, const LabelPropagationParameters& params) {
  str << "  Label Propagation Parameters:" << std::endl;
  str << "    Algorithm:                        " << params.algorithm << std::endl;
  str << "    Maximum Iterations:               " << params.maximum_iterations << std::endl;
  str << "    Part Weight Update Factor:        " << params.part_weight_update_factor << std::endl;
  str << "    Localized:                        " << std::boolalpha << params.localized << std::endl;
  str << "    Numa Aware:                       " << std::boolalpha << params.numa_aware << std::endl;
  str << "    Rebalancing:                      " << std::boolalpha << params.rebalancing << std::endl;
  if ( !params.localized ) {
    str << "    Execution Policy:                 " << params.execution_policy << std::endl;
    str << "    Execution Policy Alpha:           " << params.execution_policy_alpha << std::endl;
  }
  return str;
}

struct RefinementParameters {
  LabelPropagationParameters label_propagation;
  bool use_batch_uncontractions = false;
  size_t batch_size = 1;
};

inline std::ostream & operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  str << "  Use Batch Uncontractions:           " << std::boolalpha << params.use_batch_uncontractions << std::endl;
  if (params.use_batch_uncontractions) {
    str << "  Batch Size:                         " << params.batch_size << std::endl;
  }
  str << std::endl << params.label_propagation;
  return str;
}

struct SharedMemoryParameters {
  size_t num_threads = 1;
  size_t shuffle_block_size = 2;
  InitialHyperedgeDistribution initial_hyperedge_distribution = InitialHyperedgeDistribution::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const SharedMemoryParameters& params) {
  str << "Shared Memory Parameters:             " << std::endl;
  str << "  Number of Threads:                  " << params.num_threads << std::endl;
  str << "  Random Shuffle Block Size:          " << params.shuffle_block_size << std::endl;
  str << "  Initial Hyperedge Distribution:     " << params.initial_hyperedge_distribution << std::endl;
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
  kahypar::ContextType type = kahypar::ContextType::main;

  Context() { }

  bool isMainRecursiveBisection() const {
    return partition.mode == kahypar::Mode::recursive_bisection &&
           type == kahypar::ContextType::main;
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

  void setupContractionLimit(const HypernodeWeight total_hypergraph_weight) {
    // Setup contraction limit
    if (initial_partitioning.mode == InitialPartitioningMode::recursive) {
      coarsening.contraction_limit =
        2 * std::max(shared_memory.num_threads, static_cast<size_t>(partition.k)) *
        coarsening.contraction_limit_multiplier;
    } else {
      coarsening.contraction_limit =
        coarsening.contraction_limit_multiplier * partition.k;
    }

    // Setup maximum allowed vertex and high-degree vertex weight
    setupMaximumAllowedNodeWeight(total_hypergraph_weight);
  }

  void setupMaximumAllowedNodeWeight(const HypernodeWeight total_hypergraph_weight) {
    HypernodeWeight min_block_weight = std::numeric_limits<HypernodeWeight>::max();
    for ( PartitionID part_id = 0; part_id < partition.k; ++part_id ) {
      min_block_weight = std::min(min_block_weight, partition.max_part_weights[part_id]);
    }

    double hypernode_weight_fraction =
      coarsening.max_allowed_weight_multiplier
      / coarsening.contraction_limit;
    double high_degree_hypernode_weight_fraction =
      coarsening.max_allowed_high_degree_node_weight_multiplier
      / coarsening.contraction_limit;
    coarsening.max_allowed_node_weight =
      std::ceil(hypernode_weight_fraction * total_hypergraph_weight);
    coarsening.max_allowed_high_degree_node_weight =
      std::ceil(high_degree_hypernode_weight_fraction * total_hypergraph_weight);
    coarsening.max_allowed_node_weight =
      std::min(coarsening.max_allowed_node_weight, min_block_weight);
    coarsening.max_allowed_high_degree_node_weight =
      std::min(coarsening.max_allowed_high_degree_node_weight, min_block_weight);
  }

  void sanityCheck() {
    if (partition.objective == kahypar::Objective::cut &&
        refinement.label_propagation.algorithm == LabelPropagationAlgorithm::label_propagation_km1) {
      ALGO_SWITCH("Refinement algorithm" << refinement.label_propagation.algorithm << "only works for km1 metric."
                                         << "Do you want to use the cut version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.label_propagation.algorithm
                                         << "refiner in combination with cut metric is not possible!",
                  refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_cut);
    } else if (partition.objective == kahypar::Objective::km1 &&
               refinement.label_propagation.algorithm == LabelPropagationAlgorithm::label_propagation_cut) {
      ALGO_SWITCH("Refinement algorithm" << refinement.label_propagation.algorithm << "only works for cut metric."
                                         << "Do you want to use the km1 version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.label_propagation.algorithm
                                         << "refiner in combination with km1 metric is not possible!",
                  refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_km1);
    }

    // Sanitize based on compiler flags
    if ( coarsening.algorithm == CoarseningAlgorithm::community_coarsener ) {
      ALGO_SWITCH("Coarsening algorithm" << coarsening.algorithm << "is currently not supported."
                                          << "Do you want to switch to multilevel coarsener (Y/N)?",
                  "Coarsening with" << coarsening.algorithm
                                    << "is currently not supported!",
                  coarsening.algorithm,
                  CoarseningAlgorithm::multilevel_coarsener);
    }

    if ( !preprocessing.use_community_detection ) {
      if ( coarsening.algorithm == CoarseningAlgorithm::community_coarsener ) {
        ALGO_SWITCH("Coarsening algorithm" << coarsening.algorithm << "only works if community detection is enabled."
                                           << "Do you want to enable community detection (Y/N)?",
                    "Coarsening with" << coarsening.algorithm
                                      << "without community detection is not possible!",
                    preprocessing.use_community_detection,
                    true);
      } else if ( preprocessing.use_community_redistribution ) {
        ALGO_SWITCH("Community redistribution only works if community detection is enabled."
                    << "Do you want to enable community detection (Y/N)?",
                    "Community redistribution without community detection is not possible!",
                    preprocessing.use_community_detection,
                    true);
      }
    }

    if ( refinement.label_propagation.localized ) {
      // If we use localized label propagation, we want to execute LP on each level
      // only on the uncontracted hypernodes
      refinement.label_propagation.execution_policy = ExecutionType::constant;
      refinement.label_propagation.execution_policy_alpha = 1.0;
    }

    switch ( coarsening.algorithm ) {
      case CoarseningAlgorithm::community_coarsener:
        partition.paradigm = Paradigm::nlevel;
        break;
      case CoarseningAlgorithm::multilevel_coarsener:
        partition.paradigm = Paradigm::multilevel;
        refinement.label_propagation.execution_policy = ExecutionType::always;
        refinement.label_propagation.localized = false;
        break;
      case CoarseningAlgorithm::UNDEFINED:
        break;
    }
  }
};

inline std::ostream & operator<< (std::ostream& str, const Context& context) {
  str << "*******************************************************************************\n"
      << "*                            Partitioning Context                             *\n"
      << "*******************************************************************************\n"
      << context.partition
      << "-------------------------------------------------------------------------------\n"
      << context.preprocessing
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
}  // namespace mt_kahypar

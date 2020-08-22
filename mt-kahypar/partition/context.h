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
  #ifdef KAHYPAR_USE_N_LEVEL_PARADIGM
  Paradigm paradigm = Paradigm::nlevel;
  #else
  Paradigm paradigm = Paradigm::multilevel;
  #endif
  kahypar::Mode mode = kahypar::Mode::UNDEFINED;
  kahypar::Objective objective = kahypar::Objective::UNDEFINED;
  double epsilon = std::numeric_limits<double>::max();
  PartitionID k = std::numeric_limits<PartitionID>::max();
  int seed = 0;
  size_t num_vcycles = 0;

  int time_limit = 0;
  std::vector<HypernodeWeight> perfect_balance_part_weights;
  std::vector<HypernodeWeight> max_part_weights;
  double large_hyperedge_size_threshold_factor = std::numeric_limits<double>::max();
  HypernodeID large_hyperedge_size_threshold = std::numeric_limits<HypernodeID>::max();
  HypernodeID ignore_hyperedge_size_threshold = std::numeric_limits<HypernodeID>::max();

  bool verbose_output = true;
  bool show_detailed_timings = false;
  bool show_detailed_clustering_timings = false;
  bool measure_detailed_uncontraction_timings = false;
  bool show_memory_consumption = false;
  bool show_advanced_cut_analysis = false;
  bool enable_progress_bar = false;
  bool sp_process_output = false;
  bool csv_output = false;
  bool write_partition_file = true;

  bool enable_profiler = false;
  int snapshot_interval = std::numeric_limits<int>::max();

  std::string graph_filename { };
  std::string graph_partition_output_folder {};
  std::string graph_partition_filename { };
  std::string graph_community_filename { };
  std::string preset_file { };
};

inline std::ostream & operator<< (std::ostream& str, const PartitioningParameters& params) {
  str << "Partitioning Parameters:" << std::endl;
  str << "  Hypergraph:                         " << params.graph_filename << std::endl;
  if ( params.write_partition_file ) {
    str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
  }
  str << "  Paradigm:                           " << params.paradigm << std::endl;
  str << "  Mode:                               " << params.mode << std::endl;
  str << "  Objective:                          " << params.objective << std::endl;
  str << "  k:                                  " << params.k << std::endl;
  str << "  epsilon:                            " << params.epsilon << std::endl;
  str << "  seed:                               " << params.seed << std::endl;
  str << "  Number of V-Cycles:                 " << params.num_vcycles << std::endl;
  str << "  Ignore HE Size Threshold:           " << params.ignore_hyperedge_size_threshold << std::endl;
  str << "  Large HE Size Threshold:            " << params.large_hyperedge_size_threshold << std::endl;
  return str;
}

struct CommunityDetectionParameters {
  LouvainEdgeWeight edge_weight_function = LouvainEdgeWeight::UNDEFINED;
  uint32_t max_pass_iterations = std::numeric_limits<uint32_t>::max();
  long double min_vertex_move_fraction = std::numeric_limits<long double>::max();
  size_t vertex_degree_sampling_threshold = std::numeric_limits<size_t>::max();
};

inline std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Edge Weight Function:             " << params.edge_weight_function << std::endl;
  str << "    Maximum Louvain-Pass Iterations:  " << params.max_pass_iterations << std::endl;
  str << "    Minimum Vertex Move Fraction:     " << params.min_vertex_move_fraction << std::endl;
  str << "    Vertex Degree Sampling Threshold: " << params.vertex_degree_sampling_threshold << std::endl;
  return str;
}

struct PreprocessingParameters {
  bool stable_construction_of_incident_edges = false;
  bool use_community_detection = false;
  CommunityDetectionParameters community_detection = { };
};

inline std::ostream & operator<< (std::ostream& str, const PreprocessingParameters& params) {
  str << "Preprocessing Parameters:" << std::endl;
  str << "  Use Community Detection:            " << std::boolalpha << params.use_community_detection << std::endl;
  if (params.use_community_detection) {
    str << std::endl << params.community_detection;
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
  bool use_adaptive_edge_size = false;
  bool use_adaptive_max_allowed_node_weight = false;
  double max_allowed_weight_fraction = std::numeric_limits<double>::max();
  double adaptive_node_weight_shrink_factor_threshold = std::numeric_limits<double>::max();
  double max_allowed_weight_multiplier = std::numeric_limits<double>::max();
  double minimum_shrink_factor = std::numeric_limits<double>::max();
  double maximum_shrink_factor = std::numeric_limits<double>::max();
  size_t vertex_degree_sampling_threshold = std::numeric_limits<size_t>::max();

  // Those will be determined dynamically
  HypernodeWeight max_allowed_node_weight = 0;
  HypernodeID contraction_limit = 0;
};

inline std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params) {
  str << "Coarsening Parameters:" << std::endl;
  str << "  Algorithm:                          " << params.algorithm << std::endl;
  str << "  Use Adaptive Edge Size:             " << std::boolalpha << params.use_adaptive_edge_size << std::endl;
  #ifdef KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES
  str << "  Use Adaptive Max Node Weight:       " << std::boolalpha << params.use_adaptive_max_allowed_node_weight << std::endl;
  #endif
  if ( params.use_adaptive_max_allowed_node_weight ) {
    str << "  Max Allowed Weight Fraction:        " << params.max_allowed_weight_fraction << std::endl;
    str << "  Adaptive Node Weight Threshold:     " << params.adaptive_node_weight_shrink_factor_threshold << std::endl;
    str << "  Initial Max Hypernode Weight:       " << params.max_allowed_node_weight << std::endl;
  } else {
    str << "  Max Allowed Weight Multiplier:      " << params.max_allowed_weight_multiplier << std::endl;
    str << "  Maximum Allowed Hypernode Weight:   " << params.max_allowed_node_weight << std::endl;
  }
  str << "  Contraction Limit Multiplier:       " << params.contraction_limit_multiplier << std::endl;
  str << "  Contraction Limit:                  " << params.contraction_limit << std::endl;
  if ( params.algorithm == CoarseningAlgorithm::multilevel_coarsener ) {
    str << "  Minimum Shrink Factor:              " << params.minimum_shrink_factor << std::endl;
    str << "  Maximum Shrink Factor:              " << params.maximum_shrink_factor << std::endl;
  }
  str << "  Vertex Degree Sampling Threshold:   " << params.vertex_degree_sampling_threshold << std::endl;
  str << std::endl << params.rating;
  return str;
}

struct LabelPropagationParameters {
  LabelPropagationAlgorithm algorithm = LabelPropagationAlgorithm::do_nothing;
  size_t maximum_iterations = 1;
  bool rebalancing = true;
  bool execute_sequential = false;
  size_t hyperedge_size_activation_threshold = std::numeric_limits<size_t>::max();
};

inline std::ostream & operator<< (std::ostream& str, const LabelPropagationParameters& params) {
  str << "  Label Propagation Parameters:" << std::endl;
  str << "    Algorithm:                        " << params.algorithm << std::endl;
  if ( params.algorithm != LabelPropagationAlgorithm::do_nothing ) {
    str << "    Maximum Iterations:               " << params.maximum_iterations << std::endl;
    str << "    Rebalancing:                      " << std::boolalpha << params.rebalancing << std::endl;
    str << "    HE Size Activation Threshold:     " << std::boolalpha << params.hyperedge_size_activation_threshold << std::endl;
  }
  return str;
}

struct FMParameters {
  FMAlgorithm algorithm = FMAlgorithm::do_nothing;
  size_t multitry_rounds = 0;
  bool perform_moves_global = false;
  bool revert_parallel = true;
  double rollback_balance_violation_factor = std::numeric_limits<double>::max();
  size_t num_seed_nodes = 0;
  bool shuffle = true;
  bool obey_minimal_parallelism = false;
  double min_improvement = -1.0;
  bool release_nodes = true;
  double time_limit_factor = std::numeric_limits<double>::max();
};

inline std::ostream& operator<<(std::ostream& out, const FMParameters& params) {
  out << "  FM Parameters: \n";
  out << "    Algorithm:                        " << params.algorithm << std::endl;
  if ( params.algorithm != FMAlgorithm::do_nothing ) {
    out << "    Multitry Rounds:                  " << params.multitry_rounds << std::endl;
    out << "    Perform Moves Globally:           " << std::boolalpha << params.perform_moves_global << std::endl;
    out << "    Parallel Global Rollbacks:        " << std::boolalpha << params.revert_parallel << std::endl;
    out << "    Rollback Bal. Violation Factor:   " << params.rollback_balance_violation_factor << std::endl;
    out << "    Num Seed Nodes:                   " << params.num_seed_nodes << std::endl;
    out << "    Enable Random Shuffle:            " << std::boolalpha << params.shuffle << std::endl;
    out << "    Obey Minimal Parallelism:         " << std::boolalpha << params.obey_minimal_parallelism << std::endl;
    out << "    Minimum Improvement Factor:       " << params.min_improvement << std::endl;
    out << "    Release Nodes:                    " << std::boolalpha << params.release_nodes << std::endl;
    out << "    Time Limit Factor:                " << params.time_limit_factor << std::endl;
  }
  out << std::flush;
  return out;
}

struct RefinementParameters {
  LabelPropagationParameters label_propagation;
  FMParameters fm;
  bool refine_until_no_improvement = false;
  size_t max_batch_size = std::numeric_limits<size_t>::max();
  bool initialize_gain_cache = false;
};

inline std::ostream & operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  str << "  Refine Until No Improvement:        " << std::boolalpha << params.refine_until_no_improvement << std::endl;
  #ifdef KAHYPAR_USE_N_LEVEL_PARADIGM
  str << "  Maximum Batch Size:                 " << params.max_batch_size << std::endl;
  str << "  Initialize Gain Cache:              " << std::boolalpha << params.initialize_gain_cache << std::endl;
  #endif
  str << std::endl << params.label_propagation;
  str << "\n" << params.fm;
  return str;
}

struct SparsificationParameters {
  bool use_degree_zero_contractions = false;
  bool use_heavy_net_removal = false;
  bool use_similiar_net_removal = false;
  double hyperedge_pin_weight_fraction = 0.0;
  size_t min_hash_footprint_size = 0;
  double jaccard_threshold = 1.0;
  SimiliarNetCombinerStrategy similiar_net_combiner_strategy = SimiliarNetCombinerStrategy::UNDEFINED;
  // Those will be determined dynamically
  HypernodeWeight max_hyperedge_pin_weight = std::numeric_limits<HypernodeWeight>::max();
};

inline std::ostream & operator<< (std::ostream& str, const SparsificationParameters& params) {
  str << "Sparsification Parameters:" << std::endl;
  str << "  Use Degree-Zero HN Contractions:    " << std::boolalpha << params.use_degree_zero_contractions << std::endl;
  str << "  Use Heavy Net Removal:              " << std::boolalpha << params.use_heavy_net_removal << std::endl;
  str << "  Use Similiar Net Removal:           " << std::boolalpha << params.use_similiar_net_removal << std::endl;
  if ( params.use_heavy_net_removal ) {
    str << "  Hyperedge Pin Weight Fraction:      " << params.hyperedge_pin_weight_fraction << std::endl;
    str << "  Maximum Hyperedge Pin Weight:       " << params.max_hyperedge_pin_weight << std::endl;
  }
  if ( params.use_similiar_net_removal ) {
    str << "  Min-Hash Footprint Size:            " << params.min_hash_footprint_size << std::endl;
    str << "  Jaccard Threshold:                  " << params.jaccard_threshold << std::endl;
    str << "  Similiar Net Combiner Strategy:     " << params.similiar_net_combiner_strategy << std::endl;
  }
  return str;
}

struct InitialPartitioningParameters {
  InitialPartitioningMode mode = InitialPartitioningMode::UNDEFINED;
  RefinementParameters refinement = { };
  size_t runs = 1;
  bool use_adaptive_ip_runs = false;
  size_t min_adaptive_ip_runs = std::numeric_limits<size_t>::max();
  bool use_adaptive_epsilon = false;
  bool perform_fm_refinement = false;
  bool perform_fm_refinement_on_each_bisection = false;
  size_t lp_maximum_iterations = 1;
  size_t lp_initial_block_size = 1;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params) {
  str << "Initial Partitioning Parameters:" << std::endl;
  str << "  Initial Partitioning Mode:          " << params.mode << std::endl;
  str << "  Number of Runs:                     " << params.runs << std::endl;
  str << "  Use Adaptive IP Runs:               " << std::boolalpha << params.use_adaptive_ip_runs << std::endl;
  if ( params.use_adaptive_ip_runs ) {
    str << "  Min Adaptive IP Runs:               " << params.min_adaptive_ip_runs << std::endl;
  }
  str << "  Use Adaptive Epsilon:               " << std::boolalpha << params.use_adaptive_epsilon << std::endl;
  str << "  Perform FM Refinement:              " << std::boolalpha << params.perform_fm_refinement << std::endl;
  if ( params.perform_fm_refinement ) {
    str << "  Perform FM Refinement On Bisection: " << std::boolalpha << params.perform_fm_refinement_on_each_bisection << std::endl;
  }
  str << "  Maximum Iterations of LP IP:        " << params.lp_maximum_iterations << std::endl;
  str << "  Initial Block Size of LP IP:        " << params.lp_initial_block_size << std::endl;
  str << "\nInitial Partitioning ";
  str << params.refinement << std::endl;
  return str;
}

struct SharedMemoryParameters {
  size_t num_threads = 1;
  bool use_localized_random_shuffle = false;
  size_t shuffle_block_size = 2;
};

inline std::ostream & operator<< (std::ostream& str, const SharedMemoryParameters& params) {
  str << "Shared Memory Parameters:             " << std::endl;
  str << "  Number of Threads:                  " << params.num_threads << std::endl;
  str << "  Number of used NUMA nodes:          " << TBBNumaArena::instance().num_used_numa_nodes() << std::endl;
  str << "  Use Localized Random Shuffle:       " << std::boolalpha << params.use_localized_random_shuffle << std::endl;
  str << "  Random Shuffle Block Size:          " << params.shuffle_block_size << std::endl;
  return str;
}

class Context {
 public:
  PartitioningParameters partition { };
  PreprocessingParameters preprocessing { };
  CoarseningParameters coarsening { };
  InitialPartitioningParameters initial_partitioning { };
  RefinementParameters refinement { };
  SparsificationParameters sparsification { };
  SharedMemoryParameters shared_memory { };
  kahypar::ContextType type = kahypar::ContextType::main;

  std::string algorithm_name = "MT-KaHyPar";

  Context() { }

  bool useSparsification() const {
    return sparsification.use_degree_zero_contractions ||
           sparsification.use_heavy_net_removal ||
           sparsification.use_similiar_net_removal;
  }

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

    setupSparsificationParameters();
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
    coarsening.max_allowed_node_weight =
      std::ceil(hypernode_weight_fraction * total_hypergraph_weight);
    coarsening.max_allowed_node_weight =
      std::min(coarsening.max_allowed_node_weight, min_block_weight);
  }

  void setupSparsificationParameters() {
    if ( sparsification.use_heavy_net_removal ) {
      HypernodeWeight max_block_weight = 0;
      for ( PartitionID block = 0; block < partition.k; ++block ) {
        max_block_weight = std::max(max_block_weight, partition.max_part_weights[block]);
      }

      sparsification.max_hyperedge_pin_weight = max_block_weight /
        sparsification.hyperedge_pin_weight_fraction;
    }
  }

  void sanityCheck() {
    if (partition.objective == kahypar::Objective::cut) {
      if ( refinement.label_propagation.algorithm == LabelPropagationAlgorithm::label_propagation_km1 ) {
        ALGO_SWITCH("Refinement algorithm" << refinement.label_propagation.algorithm << "only works for km1 metric."
                                          << "Do you want to use the cut version of the label propagation refiner (Y/N)?",
                    "Partitioning with" << refinement.label_propagation.algorithm
                                          << "refiner in combination with cut metric is not possible!",
                    refinement.label_propagation.algorithm,
                    LabelPropagationAlgorithm::label_propagation_cut);
      }

      if ( refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        ALGO_SWITCH("Refinement algorithm" << refinement.fm.algorithm << "only works for km1 metric."
                                          << "Do you want to disable FM refinement (Y/N)?",
                    "Partitioning with" << refinement.fm.algorithm
                                          << "refiner in combination with cut metric is not possible!",
                    refinement.fm.algorithm,
                    FMAlgorithm::do_nothing);
      }
    } else if (partition.objective == kahypar::Objective::km1 &&
               refinement.label_propagation.algorithm == LabelPropagationAlgorithm::label_propagation_cut) {
      ALGO_SWITCH("Refinement algorithm" << refinement.label_propagation.algorithm << "only works for cut metric."
                                         << "Do you want to use the km1 version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.label_propagation.algorithm
                                         << "refiner in combination with km1 metric is not possible!",
                  refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_km1);
    }

    if (partition.objective == kahypar::Objective::cut) {
      if ( initial_partitioning.refinement.label_propagation.algorithm ==
           LabelPropagationAlgorithm::label_propagation_km1 ) {
        ALGO_SWITCH("Initial Partitioning Refinement algorithm"
                                          << initial_partitioning.refinement.label_propagation.algorithm
                                          << "only works for km1 metric."
                                          << "Do you want to use the cut version of the label propagation refiner (Y/N)?",
                    "Partitioning with" << initial_partitioning.refinement.label_propagation.algorithm
                                          << "refiner in combination with cut metric is not possible!",
                    initial_partitioning.refinement.label_propagation.algorithm,
                    LabelPropagationAlgorithm::label_propagation_cut);
      }

      if ( initial_partitioning.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        ALGO_SWITCH("Initial Partitioning Refinement algorithm"
                                          << initial_partitioning.refinement.fm.algorithm
                                          << "only works for km1 metric."
                                          << "Do you want to disable FM refinement (Y/N)?",
                    "Partitioning with" << initial_partitioning.refinement.fm.algorithm
                                          << "refiner in combination with cut metric is not possible!",
                    initial_partitioning.refinement.fm.algorithm,
                    FMAlgorithm::do_nothing);
      }
    } else if (partition.objective == kahypar::Objective::km1 &&
               initial_partitioning.refinement.label_propagation.algorithm ==
               LabelPropagationAlgorithm::label_propagation_cut) {
      ALGO_SWITCH("Initial Partitioning Refinement algorithm"
                                         << initial_partitioning.refinement.label_propagation.algorithm
                                         << "only works for cut metric."
                                         << "Do you want to use the km1 version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << initial_partitioning.refinement.label_propagation.algorithm
                                         << "refiner in combination with km1 metric is not possible!",
                  initial_partitioning.refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_km1);
    }

    if ( partition.mode == kahypar::Mode::recursive_bisection ) {
      ALGO_SWITCH("Recursive bisection mode is currently not supported."
                                         << "Do you want to use the direct k-way mode instead (Y/N)?",
                  "Recursive bisection mode is currently not supported!",
                  partition.mode,
                  kahypar::Mode::direct_kway);
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
      #ifdef KAHYPAR_ENABLE_EXPERIMENTAL_FEATURES
      << context.sparsification
      << "-------------------------------------------------------------------------------\n"
      #endif
      << context.shared_memory
      << "-------------------------------------------------------------------------------";
  return str;
}
}  // namespace mt_kahypar

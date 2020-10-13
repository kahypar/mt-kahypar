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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"

#include "kahypar/partition/context_enum_classes.h"

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

  int snapshot_interval = std::numeric_limits<int>::max();

  std::string graph_filename { };
  std::string graph_partition_output_folder {};
  std::string graph_partition_filename { };
  std::string graph_community_filename { };
  std::string preset_file { };
};

std::ostream & operator<< (std::ostream& str, const PartitioningParameters& params);

struct CommunityDetectionParameters {
  LouvainEdgeWeight edge_weight_function = LouvainEdgeWeight::UNDEFINED;
  uint32_t max_pass_iterations = std::numeric_limits<uint32_t>::max();
  long double min_vertex_move_fraction = std::numeric_limits<long double>::max();
  size_t vertex_degree_sampling_threshold = std::numeric_limits<size_t>::max();
};

std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params);

struct PreprocessingParameters {
  bool stable_construction_of_incident_edges = false;
  bool use_community_detection = false;
  CommunityDetectionParameters community_detection = { };
};

std::ostream & operator<< (std::ostream& str, const PreprocessingParameters& params);

struct RatingParameters {
  RatingFunction rating_function = RatingFunction::UNDEFINED;
  HeavyNodePenaltyPolicy heavy_node_penalty_policy = HeavyNodePenaltyPolicy::UNDEFINED;
  AcceptancePolicy acceptance_policy = AcceptancePolicy::UNDEFINED;
};

std::ostream & operator<< (std::ostream& str, const RatingParameters& params);

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


std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params);

struct LabelPropagationParameters {
  LabelPropagationAlgorithm algorithm = LabelPropagationAlgorithm::do_nothing;
  size_t maximum_iterations = 1;
  bool rebalancing = true;
  bool execute_sequential = false;
  size_t hyperedge_size_activation_threshold = std::numeric_limits<size_t>::max();
};

std::ostream & operator<< (std::ostream& str, const LabelPropagationParameters& params);

struct FMParameters {
  FMAlgorithm algorithm = FMAlgorithm::do_nothing;
  size_t multitry_rounds = 1;
  bool perform_moves_global = false;
  bool revert_parallel = true;
  double rollback_balance_violation_factor = std::numeric_limits<double>::max();
  mutable size_t num_seed_nodes = 1;
  bool shuffle = true;
  mutable bool obey_minimal_parallelism = false;
  double min_improvement = -1.0;
  bool release_nodes = true;
  double time_limit_factor = std::numeric_limits<double>::max();
};

std::ostream& operator<<(std::ostream& out, const FMParameters& params);

struct NLevelGlobalFMParameters {
  bool use_global_fm = false;   // TODO this should be renamed to something more appropriate: e.g. log_level_fm or refine_after_coarsening_pass
  bool refine_until_no_improvement = false;
  size_t num_seed_nodes = 0;
  bool obey_minimal_parallelism = false;
};

std::ostream& operator<<(std::ostream& out, const NLevelGlobalFMParameters& params);

struct RefinementParameters {
  LabelPropagationParameters label_propagation;
  FMParameters fm;
  NLevelGlobalFMParameters global_fm;
  bool refine_until_no_improvement = false;
  size_t max_batch_size = std::numeric_limits<size_t>::max();
  size_t min_border_vertices_per_thread = 0;
};

std::ostream & operator<< (std::ostream& str, const RefinementParameters& params);

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

std::ostream & operator<< (std::ostream& str, const SparsificationParameters& params);

struct InitialPartitioningParameters {
  InitialPartitioningParameters() :
    // Enable all initial partitioner per default
    enabled_ip_algos(static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED), true) { }

  InitialPartitioningMode mode = InitialPartitioningMode::UNDEFINED;
  RefinementParameters refinement = { };
  std::vector<bool> enabled_ip_algos;
  size_t runs = 1;
  bool use_adaptive_ip_runs = false;
  size_t min_adaptive_ip_runs = std::numeric_limits<size_t>::max();
  bool use_adaptive_epsilon = false;
  bool perform_refinement_on_best_partitions = false;
  size_t fm_refinment_rounds = 1;
  bool remove_degree_zero_hns_before_ip = false;
  size_t lp_maximum_iterations = 1;
  size_t lp_initial_block_size = 1;
};

std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params);

struct SharedMemoryParameters {
  size_t num_threads = 1;
  bool use_localized_random_shuffle = false;
  size_t shuffle_block_size = 2;
  double degree_of_parallelism = 1.0;
};

std::ostream & operator<< (std::ostream& str, const SharedMemoryParameters& params);

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

  std::string algorithm_name = "Mt-KaHyPar";

  Context() { }

  bool useSparsification() const ;

  bool isMainRecursiveBisection() const ;

  void setupPartWeights(const HypernodeWeight total_hypergraph_weight);

  void setupContractionLimit(const HypernodeWeight total_hypergraph_weight);

  void setupMaximumAllowedNodeWeight(const HypernodeWeight total_hypergraph_weight);

  void setupSparsificationParameters();

  void sanityCheck();
};

std::ostream & operator<< (std::ostream& str, const Context& context);

}  // namespace mt_kahypar

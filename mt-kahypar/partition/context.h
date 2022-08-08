/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/context_enum_classes.h"

#include "kahypar/partition/context_enum_classes.h"

namespace mt_kahypar {
struct PartitioningParameters {
  #ifdef USE_STRONG_PARTITIONER
  Paradigm paradigm = Paradigm::nlevel;
  #else
  Paradigm paradigm = Paradigm::multilevel;
  #endif
  Mode mode = Mode::UNDEFINED;
  kahypar::Objective objective = kahypar::Objective::UNDEFINED;
  FileFormat file_format = FileFormat::hMetis;
  InstanceType instance_type = InstanceType::UNDEFINED;
  PresetType preset_type = PresetType::UNDEFINED;
  double epsilon = std::numeric_limits<double>::max();
  PartitionID k = std::numeric_limits<PartitionID>::max();
  int seed = 0;
  size_t num_vcycles = 0;

  int time_limit = 0;
  bool use_individual_part_weights = false;
  std::vector<HypernodeWeight> perfect_balance_part_weights;
  std::vector<HypernodeWeight> max_part_weights;
  double large_hyperedge_size_threshold_factor = std::numeric_limits<double>::max();
  HypernodeID large_hyperedge_size_threshold = std::numeric_limits<HypernodeID>::max();
  HypernodeID ignore_hyperedge_size_threshold = std::numeric_limits<HypernodeID>::max();

  bool verbose_output = true;
  bool show_detailed_timings = false;
  bool show_detailed_clustering_timings = false;
  bool measure_detailed_uncontraction_timings = false;
  size_t timings_output_depth = std::numeric_limits<size_t>::max();
  bool show_memory_consumption = false;
  bool show_advanced_cut_analysis = false;
  bool enable_progress_bar = false;
  bool sp_process_output = false;
  bool csv_output = false;
  bool write_partition_file = false;
  bool deterministic = false;

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
  bool low_memory_contraction = false;
  long double min_vertex_move_fraction = std::numeric_limits<long double>::max();
  size_t vertex_degree_sampling_threshold = std::numeric_limits<size_t>::max();
  size_t num_sub_rounds_deterministic = 16;
};

std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params);

struct PreprocessingParameters {
  bool stable_construction_of_incident_edges = false;
  bool use_community_detection = false;
  bool disable_community_detection_for_mesh_graphs = true;
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
  size_t num_sub_rounds_deterministic = 16;

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
  mutable size_t num_seed_nodes = 1;

  double rollback_balance_violation_factor = std::numeric_limits<double>::max();
  double min_improvement = -1.0;
  double time_limit_factor = std::numeric_limits<double>::max();

  bool perform_moves_global = false;
  bool rollback_parallel = true;
  bool iter_moves_on_recalc = false;
  bool shuffle = true;
  mutable bool obey_minimal_parallelism = false;
  bool release_nodes = true;
};

std::ostream& operator<<(std::ostream& out, const FMParameters& params);

struct NLevelGlobalFMParameters {
  bool use_global_fm = false;   // TODO this should be renamed to something more appropriate: e.g. log_level_fm or refine_after_coarsening_pass
  bool refine_until_no_improvement = false;
  size_t num_seed_nodes = 0;
  bool obey_minimal_parallelism = false;
};

std::ostream& operator<<(std::ostream& out, const NLevelGlobalFMParameters& params);

struct FlowParameters {
  FlowAlgorithm algorithm = FlowAlgorithm::do_nothing;
  double alpha = 0.0;
  HypernodeID max_num_pins = std::numeric_limits<HypernodeID>::max();
  bool find_most_balanced_cut = false;
  bool determine_distance_from_cut = false;
  double parallel_searches_multiplier = 1.0;
  size_t num_parallel_searches = 0;
  size_t max_bfs_distance = 0;
  double min_relative_improvement_per_round = 0.0;
  double time_limit_factor = 0.0;
  bool skip_small_cuts = false;
  bool skip_unpromising_blocks = false;
  bool pierce_in_bulk = false;
};

std::ostream& operator<<(std::ostream& out, const FlowParameters& params);

struct DeterministicRefinementParameters {
  size_t num_sub_rounds_sync_lp = 5;
  bool use_active_node_set = false;
  bool recalculate_gains_on_second_apply = false;
};

std::ostream& operator<<(std::ostream& out, const DeterministicRefinementParameters& params);

struct RefinementParameters {
  LabelPropagationParameters label_propagation;
  FMParameters fm;
  DeterministicRefinementParameters deterministic_refinement;
  NLevelGlobalFMParameters global_fm;
  FlowParameters flows;
  bool refine_until_no_improvement = false;
  double relative_improvement_threshold = 0.0;
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

  Mode mode = Mode::UNDEFINED;
  RefinementParameters refinement = { };
  std::vector<bool> enabled_ip_algos;
  size_t runs = 1;
  bool use_adaptive_ip_runs = false;
  size_t min_adaptive_ip_runs = std::numeric_limits<size_t>::max();
  bool perform_refinement_on_best_partitions = false;
  size_t fm_refinment_rounds = 1;
  bool remove_degree_zero_hns_before_ip = false;
  size_t lp_maximum_iterations = 1;
  size_t lp_initial_block_size = 1;
  size_t population_size = 16;
};

std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params);

struct SharedMemoryParameters {
  size_t num_threads = 1;
  size_t static_balancing_work_packages = 128;
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
  mutable size_t initial_km1 = std::numeric_limits<size_t>::max();

  Context() { }

  bool useSparsification() const ;

  void setupPartWeights(const HypernodeWeight total_hypergraph_weight);

  void setupContractionLimit(const HypernodeWeight total_hypergraph_weight);

  void setupMaximumAllowedNodeWeight(const HypernodeWeight total_hypergraph_weight);

  void setupSparsificationParameters();

  void setupThreadsPerFlowSearch();

  void sanityCheck();
};

std::ostream & operator<< (std::ostream& str, const Context& context);

}  // namespace mt_kahypar

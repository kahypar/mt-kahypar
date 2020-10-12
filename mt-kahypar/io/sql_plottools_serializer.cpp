/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "sql_plottools_serializer.h"

#include <sstream>


#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::io::serializer {

std::string serialize(const PartitionedHypergraph& hypergraph,
                                    const Context& context,
                                    const std::chrono::duration<double>& elapsed_seconds) {
  if (context.partition.sp_process_output) {
    std::stringstream oss;
    oss << "RESULT"
        << " algorithm=" << context.algorithm_name
        << " graph=" << context.partition.graph_filename.substr(
            context.partition.graph_filename.find_last_of('/') + 1)
        << " numHNs=" << hypergraph.initialNumNodes()
        << " numHEs=" << hypergraph.initialNumEdges()
        << " paradigm=" << context.partition.paradigm
        << " mode=" << context.partition.mode
        << " objective=" << context.partition.objective
        << " k=" << context.partition.k
        << " epsilon=" << context.partition.epsilon
        << " seed=" << context.partition.seed
        << " num_vcycles=" << context.partition.num_vcycles
        << " large_hyperedge_size_threshold_factor=" << context.partition.large_hyperedge_size_threshold_factor
        << " large_hyperedge_size_threshold=" << context.partition.large_hyperedge_size_threshold
        << " ignore_hyperedge_size_threshold=" << context.partition.ignore_hyperedge_size_threshold
        << " time_limit=" << context.partition.time_limit
        << " perfect_balance_part_weight=" << context.partition.perfect_balance_part_weights[0]
        << " max_part_weight=" << context.partition.max_part_weights[0]
        << " total_graph_weight=" << hypergraph.totalWeight()
        << " use_community_detection=" << std::boolalpha << context.preprocessing.use_community_detection
        << " community_edge_weight_function=" << context.preprocessing.community_detection.edge_weight_function
        << " community_max_pass_iterations=" << context.preprocessing.community_detection.max_pass_iterations
        << " community_min_vertex_move_fraction=" << context.preprocessing.community_detection.min_vertex_move_fraction
        << " community_vertex_degree_sampling_threshold=" << context.preprocessing.community_detection.vertex_degree_sampling_threshold
        << " coarsening_algorithm=" << context.coarsening.algorithm
        << " coarsening_contraction_limit_multiplier=" << context.coarsening.contraction_limit_multiplier
        << " coarsening_use_adaptive_edge_size=" << std::boolalpha << context.coarsening.use_adaptive_edge_size
        << " coarsening_use_adaptive_max_allowed_node_weight=" << std::boolalpha << context.coarsening.use_adaptive_max_allowed_node_weight
        << " coarsening_max_allowed_weight_fraction=" << context.coarsening.max_allowed_weight_fraction
        << " coarsening_adaptive_node_weight_shrink_factor_threshold=" << context.coarsening.adaptive_node_weight_shrink_factor_threshold
        << " coarsening_max_allowed_weight_multiplier=" << context.coarsening.max_allowed_weight_multiplier
        << " coarsening_minimum_shrink_factor=" << context.coarsening.minimum_shrink_factor
        << " coarsening_maximum_shrink_factor=" << context.coarsening.maximum_shrink_factor
        << " coarsening_max_allowed_node_weight=" << context.coarsening.max_allowed_node_weight
        << " coarsening_vertex_degree_sampling_threshold=" << context.coarsening.vertex_degree_sampling_threshold
        << " coarsening_contraction_limit=" << context.coarsening.contraction_limit
        << " rating_function=" << context.coarsening.rating.rating_function
        << " rating_heavy_node_penalty_policy=" << context.coarsening.rating.heavy_node_penalty_policy
        << " rating_acceptance_policy=" << context.coarsening.rating.acceptance_policy
        << " initial_partitioning_mode=" << context.initial_partitioning.mode
        << " initial_partitioning_runs=" << context.initial_partitioning.runs
        << " initial_partitioning_use_adaptive_ip_runs=" << std::boolalpha << context.initial_partitioning.use_adaptive_ip_runs
        << " initial_partitioning_min_adaptive_ip_runs=" << context.initial_partitioning.min_adaptive_ip_runs
        << " initial_partitioning_use_adaptive_epsilon=" << std::boolalpha << context.initial_partitioning.use_adaptive_epsilon
        << " initial_partitioning_perform_refinement_on_best_partitions=" << std::boolalpha << context.initial_partitioning.perform_refinement_on_best_partitions
        << " initial_partitioning_fm_refinment_rounds=" << std::boolalpha << context.initial_partitioning.fm_refinment_rounds
        << " initial_partitioning_remove_degree_zero_hns_before_ip=" << std::boolalpha << context.initial_partitioning.remove_degree_zero_hns_before_ip
        << " initial_partitioning_lp_maximum_iterations=" << context.initial_partitioning.lp_maximum_iterations
        << " initial_partitioning_lp_initial_block_size=" << context.initial_partitioning.lp_initial_block_size
        << " sparsification_use_degree_zero_contractions=" << std::boolalpha << context.sparsification.use_degree_zero_contractions
        << " sparsification_use_heavy_net_removal=" << std::boolalpha << context.sparsification.use_heavy_net_removal
        << " sparsification_use_similiar_net_removal=" << std::boolalpha << context.sparsification.use_similiar_net_removal
        << " sparsification_hyperedge_pin_weight_fraction=" << context.sparsification.hyperedge_pin_weight_fraction
        << " sparsification_max_hyperedge_pin_weight=" << context.sparsification.max_hyperedge_pin_weight
        << " sparsification_min_hash_footprint_size=" << context.sparsification.min_hash_footprint_size
        << " sparsification_jaccard_threshold=" << context.sparsification.jaccard_threshold
        << " sparsification_similiar_net_combiner_strategy=" << context.sparsification.similiar_net_combiner_strategy
        << " refine_until_no_improvement=" << std::boolalpha << context.refinement.refine_until_no_improvement
        << " max_batch_size=" << context.refinement.max_batch_size
        << " min_border_vertices_per_thread=" << context.refinement.min_border_vertices_per_thread
        << " lp_algorithm=" << context.refinement.label_propagation.algorithm
        << " lp_maximum_iterations=" << context.refinement.label_propagation.maximum_iterations
        << " lp_rebalancing=" << std::boolalpha << context.refinement.label_propagation.rebalancing
        << " lp_hyperedge_size_activation_threshold=" << context.refinement.label_propagation.hyperedge_size_activation_threshold
        << " fm_algorithm=" << context.refinement.fm.algorithm
        << " fm_multitry_rounds=" << context.refinement.fm.multitry_rounds
        << " fm_perform_moves_global=" << std::boolalpha << context.refinement.fm.perform_moves_global
        << " fm_revert_parallel=" << std::boolalpha << context.refinement.fm.revert_parallel
        << " fm_rollback_balance_violation_factor=" << context.refinement.fm.rollback_balance_violation_factor
        << " fm_min_improvement=" << context.refinement.fm.min_improvement
        << " fm_release_nodes=" << context.refinement.fm.release_nodes
        << " fm_num_seed_nodes=" << context.refinement.fm.num_seed_nodes
        << " fm_time_limit_factor=" << context.refinement.fm.time_limit_factor
        << " fm_obey_minimal_parallelism=" << std::boolalpha << context.refinement.fm.obey_minimal_parallelism
        << " fm_shuffle=" << std::boolalpha << context.refinement.fm.shuffle
        << " global_fm_use_global_fm=" << std::boolalpha << context.refinement.global_fm.use_global_fm
        << " global_fm_refine_until_no_improvement=" << std::boolalpha << context.refinement.global_fm.refine_until_no_improvement
        << " global_fm_num_seed_nodes=" << context.refinement.global_fm.num_seed_nodes
        << " global_fm_obey_minimal_parallelism=" << std::boolalpha << context.refinement.global_fm.obey_minimal_parallelism
        << " num_threads=" << context.shared_memory.num_threads
        << " use_localized_random_shuffle=" << std::boolalpha << context.shared_memory.use_localized_random_shuffle
        << " shuffle_block_size=" << context.shared_memory.shuffle_block_size;

    // Metrics
    if ( hypergraph.initialNumEdges() > 0 ) {
      oss << " cut=" << metrics::hyperedgeCut(hypergraph)
          << " soed=" << metrics::soed(hypergraph)
          << " km1=" << metrics::km1(hypergraph)
          << " imbalance=" << metrics::imbalance(hypergraph, context);
    }
    oss << " totalPartitionTime=" << elapsed_seconds.count();

    // Timings
    utils::Timer::instance(context.partition.show_detailed_timings).serialize(oss);

    // Stats
    oss << utils::Stats::instance();

    // Initial Partitioning Stats
    oss << utils::InitialPartitioningStats::instance();

    return oss.str();
  } else {
    return "";
  }
}
}  // namespace

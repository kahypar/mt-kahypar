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

#pragma once
#include <array>
#include <chrono>
#include <sstream>
#include <string>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace io {
namespace serializer {
static inline std::string serialize(const Hypergraph& hypergraph,
                                    const Context& context,
                                    const std::chrono::duration<double>& elapsed_seconds) {
  if (context.partition.sp_process_output) {
    std::stringstream oss;
    oss << "RESULT"
        << " graph=" << context.partition.graph_filename.substr(
      context.partition.graph_filename.find_last_of('/') + 1)
        << " numHNs=" << hypergraph.initialNumNodes()
        << " numHEs=" << hypergraph.initialNumEdges()
        << " mode=" << context.partition.mode
        << " objective=" << context.partition.objective
        << " k=" << context.partition.k
        << " epsilon=" << context.partition.epsilon
        << " seed=" << context.partition.seed
        << " hyperedge_size_threshold=" << context.partition.hyperedge_size_threshold
        << " time_limit=" << context.partition.time_limit
        << " perfect_balance_part_weight=" << context.partition.perfect_balance_part_weights[0]
        << " max_part_weight=" << context.partition.max_part_weights[0]
        << " total_graph_weight=" << hypergraph.totalWeight()
        << " use_community_structure_from_file=" << std::boolalpha << context.preprocessing.use_community_structure_from_file
        << " community_load_balancing_strategy=" << context.preprocessing.community_detection.load_balancing_strategy
        << " community_size_constraint_factor=" << context.preprocessing.community_detection.size_constraint_factor
        << " community_edge_weight_function=" << context.preprocessing.community_detection.edge_weight_function
        << " community_max_pass_iterations=" << context.preprocessing.community_detection.max_pass_iterations
        << " community_min_eps_improvement=" << context.preprocessing.community_detection.min_eps_improvement
        << " use_community_redistribution=" << std::boolalpha << context.preprocessing.community_redistribution.use_community_redistribution
        << " community_redistribution_assignment_strategy=" << context.preprocessing.community_redistribution.assignment_strategy
        << " community_redistribution_assignment_objective=" << context.preprocessing.community_redistribution.assignment_objective
        << " coarsening_algorithm=" << context.coarsening.algorithm
        << " coarsening_contraction_limit_multiplier=" << context.coarsening.contraction_limit_multiplier
        << " coarsening_max_allowed_weight_multiplier=" << context.coarsening.max_allowed_weight_multiplier
        << " coarsening_max_allowed_high_degree_node_weight_multiplier=" << context.coarsening.max_allowed_high_degree_node_weight_multiplier
        << " coarsening_use_high_degree_vertex_threshold=" << std::boolalpha << context.coarsening.use_high_degree_vertex_threshold
        << " coarsening_max_allowed_node_weight=" << context.coarsening.max_allowed_node_weight
        << " coarsening_max_allowed_high_degree_node_weight=" << context.coarsening.max_allowed_high_degree_node_weight
        << " coarsening_contraction_limit=" << context.coarsening.contraction_limit
        << " coarsening_high_degree_vertex_threshold=" << context.coarsening.high_degree_vertex_threshold
        << " rating_function=" << context.coarsening.rating.rating_function
        << " rating_heavy_node_penalty_policy=" << context.coarsening.rating.heavy_node_penalty_policy
        << " rating_acceptance_policy=" << context.coarsening.rating.acceptance_policy
        << " initial_partitioning_mode=" << context.initial_partitioning.mode
        << " initial_partitioning_runs=" << context.initial_partitioning.runs
        << " initial_partitioning_use_adaptive_epsilon=" << std::boolalpha << context.initial_partitioning.use_adaptive_epsilon
        << " initial_partitioning_lp_maximum_iterations=" << context.initial_partitioning.lp_maximum_iterations
        << " initial_partitioning_lp_initial_block_size=" << context.initial_partitioning.lp_initial_block_size
        << " use_batch_uncontractions=" << std::boolalpha << context.refinement.use_batch_uncontractions
        << " batch_size=" << context.refinement.batch_size
        << " lp_algorithm=" << context.refinement.label_propagation.algorithm
        << " lp_maximum_iterations=" << context.refinement.label_propagation.maximum_iterations
        << " lp_part_weight_update_factor=" << context.refinement.label_propagation.part_weight_update_factor
        << " lp_localized=" << std::boolalpha << context.refinement.label_propagation.localized
        << " lp_numa_aware=" << std::boolalpha << context.refinement.label_propagation.numa_aware
        << " lp_rebalancing=" << std::boolalpha << context.refinement.label_propagation.rebalancing
        << " lp_execution_policy=" << context.refinement.label_propagation.execution_policy
        << " lp_execution_policy_alpha=" << context.refinement.label_propagation.execution_policy_alpha
        << " num_threads=" << context.shared_memory.num_threads
        << " initial_hyperedge_distribution=" << context.shared_memory.initial_hyperedge_distribution;

    // Metrics
    if ( hypergraph.initialNumEdges() > 0 ) {
      oss << " cut=" << metrics::hyperedgeCut(hypergraph)
          << " soed=" << metrics::soed(hypergraph)
          << " km1=" << metrics::km1(hypergraph)
          << " absorption=" << metrics::absorption(hypergraph)
          << " imbalance=" << metrics::imbalance(hypergraph, context);
    }
    oss << " totalPartitionTime=" << elapsed_seconds.count();

    // Part Weights and Sizes
    /*for (PartitionID i = 0; i < context.partition.k; ++i) {
      oss << " partSize" << i << "=" << hypergraph.partSize(i);
    }
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      oss << " partWeight" << i << "=" << hypergraph.partWeight(i);
    }*/

    // Timings
    utils::Timer::instance(context.partition.detailed_timings).serialize(oss);

    // Stats
    oss << utils::Stats::instance();

    // Initial Partitioning Stats
    oss << utils::InitialPartitioningStats::instance();

    return oss.str();
  }
  return "";
}
}  // namespace serializer
}  // namespace io
}  // namespace mt_kahypar

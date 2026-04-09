/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2017 Robin Andre <robinandre1995@web.de>
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

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/evolutionary/edge_frequency.h"
#include "mt-kahypar/partition/evolutionary/population.h"



namespace mt_kahypar::combine {
static constexpr bool debug = false;

  vec<PartitionID> combinePartitions(const std::vector<std::vector<PartitionID>> &parent_partitions);

  template <typename TypeTraits>
  Individual usingKWaySelection(const typename TypeTraits::Hypergraph& input_hg, TargetGraph* target_graph, Population& population, Context context, std::mt19937* rng) {
    std::vector<size_t> parents;
    //Maybe change to actually use Tournament Selection
    size_t best(population.randomIndividualSafe(context, rng));
    parents.push_back(best);
    for (int x = 1; x < context.evolutionary.kway_combine; x++) {
      size_t new_parent = population.randomIndividualSafe(context, rng);
      parents.push_back(new_parent);
      if (population.fitnessAtSafe(new_parent) <= population.fitnessAtSafe(best)) {
        best = new_parent;
      }
    }

    std::vector<PartitionID> best_partition = population.partitionCopySafe(best);
    std::unordered_map<PartitionID, int> comm_to_block;

    // aquire lock --- possibly unnecessary
    std::vector<std::vector<PartitionID>> parent_partitions;
    for (auto parent_id : parents) {
      parent_partitions.push_back(population.partitionCopySafe(parent_id)); // FIXED
    }

    vec<PartitionID> comms = combinePartitions(parent_partitions);

    typename TypeTraits::Hypergraph hypergraph = input_hg.copy(parallel_tag_t{});
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);

    for (const HypernodeID& hn : hypergraph.nodes()) {
      partitioned_hypergraph.setOnlyNodePart(hn, best_partition[hn]);
      if (comm_to_block.find(comms[hn]) == comm_to_block.end()) {
        comm_to_block[comms[hn]] = best_partition[hn];
      }
    }

    partitioned_hypergraph.initializePartition();
    hypergraph.setCommunityIDs(std::move(comms));

    if (context.partition.mode == Mode::direct) {
      Context vc_context(context);
      vc_context.setupPartWeights(hypergraph.totalWeight());
      Multilevel<TypeTraits>::evolutionPartitionVCycle(
          hypergraph, partitioned_hypergraph, vc_context, comm_to_block, target_graph);
    } else {
      throw InvalidParameterException("Invalid partitioning mode!");
    }

    return Individual(partitioned_hypergraph, context);



    return Individual(2);
  }

/*template<typename Hypergraph>
Individual usingTournamentSelection(Hypergraph& hg, const Context& context, const Population& population) {
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>() };
  temporary_context.coarsening.rating.rating_function = RatingFunction::heavy_edge;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::evolutionary;

  const auto& parents = population.tournamentSelect();


  return combine::partitions(hg, parents, temporary_context);
}
template<typename Hypergraph>
Individual edgeFrequency(Hypergraph& hg, const Context& context, const Population& population) {
  const HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  hg.reset();
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>(),
             meta::Int2Type<static_cast<int>(EvoCombineStrategy::edge_frequency)>() };

  temporary_context.coarsening.rating.rating_function = RatingFunction::edge_frequency;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::normal;
  temporary_context.coarsening.rating.heavy_node_penalty_policy =
    HeavyNodePenaltyPolicy::edge_frequency_penalty;

  temporary_context.evolutionary.edge_frequency =
    computeEdgeFrequency(population.listOfBest(context.evolutionary.edge_frequency_amount),
                         hg.initialNumEdges());

  DBG << V(temporary_context.evolutionary.action.decision());


  Partitioner().partition(hg, temporary_context);

  const HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::evolutionary,
                        std::chrono::duration<double>(end - start).count());


  DBG << "final result" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
  io::serializer::serializeEvolutionary(temporary_context, hg);
  return Individual(hg, context);
}*/
} // namespace mt_kahypar::combine


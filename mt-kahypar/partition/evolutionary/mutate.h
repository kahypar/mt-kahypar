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

#include <vector>

#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/evolutionary/edge_frequency.h"
#include "mt-kahypar/partition/partitioner.h"

namespace mt_kahypar::mutate {
  static constexpr bool debug = false;

  template<typename TypeTraits>
  Individual vCycle(typename TypeTraits::Hypergraph& hypergraph, std::vector<PartitionID> cur, TargetGraph* target_graph, Context context) {
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);

    vec<PartitionID> comms(hypergraph.initialNumNodes());
    std::unordered_map<PartitionID, int> comm_to_block;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      partitioned_hypergraph.setOnlyNodePart(hn, cur[hn]);
      comms[hn] = cur[hn];
    }
    for (PartitionID i = 0; i < context.partition.k; i++) {
      comm_to_block[i] = i;
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
  }

  template<typename TypeTraits>
  Individual vCycleWithNewInitialPartitioning(typename TypeTraits::Hypergraph& hypergraph, std::vector<PartitionID> cur, TargetGraph* target_graph, Context context) {
    vec<PartitionID> comms(hypergraph.initialNumNodes());
    for (const HypernodeID& hn : hypergraph.nodes()) {
      comms[hn] = cur[hn];
    }
    hypergraph.setCommunityIDs(std::move(comms));
    Context mut_context(context);
    if (!mut_context.partition.use_individual_part_weights) {
      mut_context.partition.max_part_weights.clear();
    }
    typename TypeTraits::PartitionedHypergraph partitioned_hypergraph = Partitioner<TypeTraits>::partition(
        hypergraph, mut_context, target_graph);

    return Individual(partitioned_hypergraph, context);
  }
}

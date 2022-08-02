/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "star_partitioning.h"

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/star_partitioning/approximate.h"
#include "mt-kahypar/partition/star_partitioning/simple_greedy.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;
using ds::SeparatedNodes;

HyperedgeWeight partition(PartitionedHypergraph& hypergraph, const Context& context, bool parallel) {
  ASSERT([&]() {
      for (const HypernodeID& hn : hypergraph.nodes()) {
        if (hypergraph.partID(hn) == kInvalidPartition) {
          return false;
        }
      }
      return true;
    } (), "There are unassigned hypernodes!");

  Array<HypernodeWeight> part_weights(hypergraph.k());
  for (PartitionID part = 0; part < hypergraph.k(); ++part) {
    part_weights[part] = hypergraph.partWeight(part);
  }

  if (context.partition.star_partitioning_algorithm == StarPartitioningAlgorithm::simple_greedy) {
    SimpleGreedy sg(context.partition.k);
    sg.partition(hypergraph, context, part_weights, parallel);
  } else if (context.partition.star_partitioning_algorithm == StarPartitioningAlgorithm::approximate) {
    Approximate ap(context.partition.k);
    if (parallel) {
      ap.partition(hypergraph, context, part_weights, parallel_tag_t());
    } else {
      ap.partition(hypergraph, context, part_weights);
    }
  }

  const SeparatedNodes& sn = hypergraph.separatedNodes();
  if (parallel) {
    tbb::enumerable_thread_specific<HyperedgeWeight> cut(0);
    tbb::parallel_for(ID(0), sn.numNodes(), [&](const HypernodeID node) {
      HyperedgeWeight& local_sum = cut.local();
      for (const auto& e: sn.inwardEdges(node)) {
        if (hypergraph.partID(e.target) != hypergraph.separatedPartID(node)) {
          local_sum += e.weight;
        }
      }
    });
    return cut.combine(std::plus<>());
  } else {
    HyperedgeWeight cut = 0;
    for (HypernodeID node = 0; node < sn.numNodes(); ++node) {
      for (const auto& e: sn.inwardEdges(node)) {
        if (hypergraph.partID(e.target) != hypergraph.separatedPartID(node)) {
          cut += e.weight;
        }
      }
    }
    return cut;
  }
}

void getEdgeWeightsOfNode(PartitionedHypergraph& phg, Array<HyperedgeWeight>& weights_per_part,
                          const HypernodeID& node, const HypernodeID& index) {
  const SeparatedNodes& s_nodes = phg.separatedNodes();
  for (const auto& e: s_nodes.inwardEdges(node)) {
    const PartitionID target_part = phg.partID(e.target);
    ASSERT(!phg.nodeIsEnabled(e.target) || target_part != kInvalidPartition);
    if (target_part != kInvalidPartition) {
      weights_per_part[index + target_part] += e.weight;
    }
  }
}

} // namepace star_partitioning
} // namepace mt_kahypar

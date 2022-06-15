/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/partition/initial_partitioning/recursive_bipartitioning_initial_partitioner.h"


#include <algorithm>
#include <vector>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/multilevel.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"

#include <tbb/parallel_sort.h>

namespace mt_kahypar {
  using ds::SeparatedNodes;
  using ds::Array;

  RecursiveBipartitioningInitialPartitioner::RecursiveBipartitioningInitialPartitioner(PartitionedHypergraph& hypergraph,
                                                                             const Context& context) :
    _hg(hypergraph),
    _context(context) { }

  void RecursiveBipartitioningInitialPartitioner::initialPartitionImpl() {
    recursive_bipartitioning::partition(_hg, _context);

    if (_context.partition.separated_nodes_processing_after_ip) {
      SeparatedNodes& separated_nodes = _hg.separatedNodes();
      tbb::enumerable_thread_specific<Array<HyperedgeWeight>> edge_weights(_hg.k());
      const HypernodeID num_s_nodes = separated_nodes.numNodes();
      Array<HyperedgeWeight> max_gains(num_s_nodes);

      tbb::parallel_for(ID(0), num_s_nodes, [&](const HypernodeID s_node) {
        Array<HyperedgeWeight>& local_edge_weights = edge_weights.local();
        local_edge_weights.assign(_hg.k(), 0);

        for (const auto& e: separated_nodes.inwardEdges(s_node)) {
          const PartitionID target_part = _hg.partID(e.target);
          if (target_part != kInvalidPartition) {
            local_edge_weights[target_part] += e.weight;
          }
        }

        HyperedgeWeight max_gain = 0;
        for (PartitionID part = 0; part < _hg.k(); ++part) {
          if (local_edge_weights[part] >= max_gain) {
            max_gain = local_edge_weights[part];
          }
        }
        max_gains[s_node] = max_gain;
      });

      Array<HypernodeID> id_order(num_s_nodes);
      tbb::parallel_for(ID(0), num_s_nodes, [&](const HypernodeID s_node) {
        id_order[s_node] = s_node;
      });

      tbb::parallel_sort(id_order.begin(), id_order.end(),
        [&](const HypernodeID& left, const HypernodeID& right) {
          return max_gains[left] > max_gains[right];
        }
      );

      Array<HypernodeWeight> part_weights(_hg.k());
      for (PartitionID part = 0; part < _hg.k(); ++part) {
        part_weights[part] = _hg.partWeight(part);
      }
      for (size_t i = 0; i < id_order.size(); ++i) {
        const HypernodeID s_node = id_order[i];
        Array<HyperedgeWeight>& local_edge_weights = edge_weights.local();
        local_edge_weights.assign(_hg.k(), 0);

        for (const auto& e: separated_nodes.inwardEdges(s_node)) {
          const PartitionID target_part = _hg.partID(e.target);
          if (target_part != kInvalidPartition) {
            local_edge_weights[target_part] += e.weight;
          }
        }

        // greedily assign separated nodes
        PartitionID max_part = kInvalidPartition;
        HyperedgeWeight max_gain = 0;
        for (PartitionID part = 0; part < _hg.k(); ++part) {
          const HypernodeWeight max_part_weight = _context.partition.max_part_weights[part];
          if (local_edge_weights[part] >= max_gain &&
              part_weights[part] + separated_nodes.nodeWeight(s_node) <= max_part_weight) {
            max_part = part;
            max_gain = local_edge_weights[part];
          }
        }
        separated_nodes.setPartID(s_node, max_part);
        part_weights[max_part] += separated_nodes.nodeWeight(s_node);
      }

      _hg.updateBlockWeights();
    }
  }
} // namepace mt_kahypar
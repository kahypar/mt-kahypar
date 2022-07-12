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
#include "mt-kahypar/partition/star_partitioning/simple_greedy.h"

namespace mt_kahypar {
  using ds::SeparatedNodes;
  using ds::Array;
  using star_partitioning::SimpleGreedy;

  RecursiveBipartitioningInitialPartitioner::RecursiveBipartitioningInitialPartitioner(PartitionedHypergraph& hypergraph,
                                                                             const Context& context) :
    _hg(hypergraph),
    _context(context) { }

  void RecursiveBipartitioningInitialPartitioner::initialPartitionImpl() {
    recursive_bipartitioning::partition(_hg, _context);

    if (_context.partition.separated_nodes_processing_after_ip) {
      SimpleGreedy sg(_context);
      SeparatedNodes& separated_nodes = _hg.separatedNodes();
      Array<HypernodeWeight> part_weights(_hg.k());

      for (PartitionID part = 0; part < _hg.k(); ++part) {
        part_weights[part] = _hg.partWeight(part);
      }

      sg.partition(separated_nodes.numNodes(), part_weights, _context.partition.max_part_weights,
        [&](Array<HyperedgeWeight>& weights, const HypernodeID node) {
          for (const auto& e: separated_nodes.inwardEdges(node)) {
            const PartitionID target_part = _hg.partID(e.target);
            if (target_part != kInvalidPartition) {
              weights[target_part] += e.weight;
            }
          }
        },
        [&](const HypernodeID node) {
          return separated_nodes.nodeWeight(node);
        },
        [&](const HypernodeID node, const PartitionID part) {
          separated_nodes.setPartID(node, part);
        });

      _hg.updateBlockWeights();
    }
  }
} // namepace mt_kahypar
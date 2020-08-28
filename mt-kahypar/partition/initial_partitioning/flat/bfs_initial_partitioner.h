/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"
#include "mt-kahypar/parallel/stl/scalable_queue.h"

namespace mt_kahypar {

class BFSInitialPartitioner : public tbb::task {
  using Queue = parallel::scalable_queue<HypernodeID>;

  static constexpr bool debug = false;

 public:
  BFSInitialPartitioner(const InitialPartitioningAlgorithm,
                         InitialPartitioningDataContainer& ip_data,
                         const Context& context,
                         const int seed) :
    _ip_data(ip_data),
    _context(context),
    _rng(seed) { }

  tbb::task* execute() override ;

 private:
  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.perfect_balance_part_weights[block];
  }

  // ! Pushes all adjacent hypernodes (not visited before) of hypernode hn
  // ! into the BFS queue of the corresponding block.
  inline void pushIncidentHypernodesIntoQueue(const PartitionedHypergraph& hypergraph,
                                              const Context& context,
                                              Queue& queue,
                                              kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                              kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                              const HypernodeID hn,
                                              const PartitionID block);

  inline void markHypernodeAsInQueue(const PartitionedHypergraph& hypergraph,
                                     kahypar::ds::FastResetFlagArray<>& hypernodes_in_queue,
                                     const HypernodeID hn,
                                     const PartitionID block) {
    hypernodes_in_queue.set(block * hypergraph.initialNumNodes() + hn, true);
  }

  inline void markHyperedgeAsInQueue(const PartitionedHypergraph& hypergraph,
                                     kahypar::ds::FastResetFlagArray<>& hyperedges_in_queue,
                                     const HyperedgeID he,
                                     const PartitionID block) {
    hyperedges_in_queue.set(block * hypergraph.initialNumEdges() + he, true);
  }

  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
  std::mt19937 _rng;
};


} // namespace mt_kahypar

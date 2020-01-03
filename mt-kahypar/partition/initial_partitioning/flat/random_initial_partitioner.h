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
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
template<typename TypeTraits>
class RandomInitialPartitionerT : public tbb::task {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<TypeTraits>;

  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;

 public:
  RandomInitialPartitionerT(InitialPartitioningDataContainer& ip_data,
                            const Context& context) :
    _ip_data(ip_data),
    _context(context) { }

  tbb::task* execute() override {
    HyperGraph& hg = _ip_data.local_hypergraph();
    int cpu_id = sched_getcpu();

    for ( const HypernodeID& hn : hg.nodes() ) {
      // Randomly select a block to assign the hypernode
      PartitionID block = utils::Randomize::instance().getRandomInt(
        0, _context.partition.k - 1, cpu_id);
      PartitionID current_block = block;
      while ( !fitsIntoBlock(hg, hn, current_block) ) {
        // If the hypernode does not fit into the random selected block
        // (because it would violate the balance constraint), we try to
        // assign it to the next block.
        current_block = ( current_block + 1 ) % _context.partition.k;
        if ( current_block == block ) {
          // In case, we find no valid block to assign the current hypernode
          // to, we assign it to random selected block (Note, will violate
          // balance constraint)
          break;
        }
      }
      hg.setNodePart(hn, current_block);
    }

    _ip_data.commit(InitialPartitioningAlgorithm::random);
    return nullptr;
  }

 private:
  bool fitsIntoBlock(HyperGraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.localPartWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.max_part_weights[block];
  }

  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
};
template <typename TypeTraits>
PartitionID RandomInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;

using RandomInitialPartitioner = RandomInitialPartitionerT<GlobalTypeTraits>;

} // namespace mt_kahypar

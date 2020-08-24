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
class RandomInitialPartitioner : public tbb::task {

  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;

 public:
  RandomInitialPartitioner(const InitialPartitioningAlgorithm,
                            InitialPartitioningDataContainer& ip_data,
                            const Context& context) :
    _ip_data(ip_data),
    _context(context) { }

  tbb::task* execute() override {
    if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::random) ) {
      HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
      PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
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
            // to, we assign it to random selected block
            break;
          }
        }
        hg.setNodePart(hn, current_block);
      }

      HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration<double>(end - start).count();
      _ip_data.commit(InitialPartitioningAlgorithm::random, time);
    }
    return nullptr;
  }

 private:
  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.perfect_balance_part_weights[block];
  }

  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
};

PartitionID RandomInitialPartitioner::kInvalidPartition = -1;

} // namespace mt_kahypar

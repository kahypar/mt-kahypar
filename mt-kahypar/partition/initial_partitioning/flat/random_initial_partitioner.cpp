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

#include "mt-kahypar/partition/initial_partitioning/flat/random_initial_partitioner.h"

#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

tbb::task* RandomInitialPartitioner::execute() {
  if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::random) ) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
    std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);

    for ( const HypernodeID& hn : hg.nodes() ) {
      // Randomly select a block to assign the hypernode
      PartitionID block = select_random_block(_rng);
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

} // namespace mt_kahypar

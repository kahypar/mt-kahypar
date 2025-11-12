#include "mt-kahypar/partition/initial_partitioning/streaming_initial_partitioner.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/policies/pseudo_peripheral_start_nodes.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

template <typename TypeTraits>
void StreamingInitialPartitioner<TypeTraits>::partitionImpl() {
  if (_ip_data.should_initial_partitioner_run(
          InitialPartitioningAlgorithm::streaming)) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();


    LOG << "Computing a" << _context.partition.k << "-way partition!";

    // This is just copied from the RandomInitialPartitioner -- replace it with your own code!
    // VVVVV
    std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);

    _ip_data.preassignFixedVertices(hg);
    for ( const HypernodeID& hn : hg.nodes() ) {
      if ( !hg.isFixed(hn) ) {
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
    }
    // ^^^^^
    // End of random initial partitioning


    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::streaming, _rng, _tag, time);
  }
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(StreamingInitialPartitioner)

} // namespace mt_kahypar

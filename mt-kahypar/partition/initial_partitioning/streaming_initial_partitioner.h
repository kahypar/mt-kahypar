#pragma once

#include "mt-kahypar/parallel/stl/scalable_queue.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"

namespace mt_kahypar {

template <typename TypeTraits>
class StreamingInitialPartitioner : public IInitialPartitioner {
  static constexpr bool debug = false;

  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

public:
  StreamingInitialPartitioner(const InitialPartitioningAlgorithm,
                              ip_data_container_t *ip_data,
                              const Context &context, const int seed,
                              const int tag)
      : _ip_data(ip::to_reference<TypeTraits>(ip_data)), _context(context),
        _rng(seed), _tag(tag) {}

private:
  void partitionImpl() final;

  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.perfect_balance_part_weights[block];
  }

  InitialPartitioningDataContainer<TypeTraits> &_ip_data;
  const Context &_context;
  std::mt19937 _rng;
  const int _tag;
};

} // namespace mt_kahypar

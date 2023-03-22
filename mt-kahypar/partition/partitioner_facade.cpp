/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/partitioner_facade.h"

#include "mt-kahypar/one_definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

namespace internal {

  template<typename TypeTraits>
  mt_kahypar_partitioned_hypergraph_t partition(mt_kahypar_hypergraph_t hypergraph,
                                                Context& context) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    Hypergraph& hg = utils::cast<Hypergraph>(hypergraph);
    PartitionedHypergraph* phg = new PartitionedHypergraph(context.partition.k, hg, parallel_tag_t { });

    // Setup Context
    context.sanityCheck();
    context.setupPartWeights(hg.totalWeight());
    context.setupContractionLimit(hg.totalWeight());
    context.setupThreadsPerFlowSearch();
    io::printContext(context);

    // Compute random partition
    utils::Randomize& rand = utils::Randomize::instance();
    phg->doParallelForAllNodes([&](const HypernodeID& hn) {
      phg->setOnlyNodePart(hn, rand.getRandomInt(0, context.partition.k - 1, SCHED_GETCPU));
    });
    phg->initializePartition();

    return mt_kahypar_partitioned_hypergraph_t {
      reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(phg), PartitionedHypergraph::TYPE };
  }

  mt_kahypar_partition_type_t getPartitionedHypergraphType(const mt_kahypar_hypergraph_t hypergraph,
                                                           const Context& context) {
    if ( hypergraph.type == STATIC_GRAPH ) {
      return STATIC_PARTITIONED_GRAPH;
    } else if ( hypergraph.type == DYNAMIC_GRAPH ) {
      return DYNAMIC_PARTITIONED_GRAPH;
    } else if ( hypergraph.type == STATIC_HYPERGRAPH ) {
      if ( context.partition.preset_type == PresetType::large_k ) {
        return STATIC_SPARSE_PARTITIONED_HYPERGRAPH;
      } else {
        return STATIC_PARTITIONED_HYPERGRAPH;
      }
    } else if ( hypergraph.type == DYNAMIC_HYPERGRAPH ) {
      return DYNAMIC_PARTITIONED_HYPERGRAPH;
    }
    return NULLPTR_PARTITION;
  }

} // namespace internal

  mt_kahypar_partitioned_hypergraph_t PartitionerFacade::partition(mt_kahypar_hypergraph_t hypergraph,
                                                                   Context& context) {
    const mt_kahypar_partition_type_t type = internal::getPartitionedHypergraphType(hypergraph, context);
    switch ( type ) {
      case STATIC_PARTITIONED_GRAPH:
        return internal::partition<StaticGraphTypeTraits>(hypergraph, context);
      case DYNAMIC_PARTITIONED_GRAPH:
        return internal::partition<DynamicGraphTypeTraits>(hypergraph, context);
      case STATIC_PARTITIONED_HYPERGRAPH:
        return internal::partition<StaticHypergraphTypeTraits>(hypergraph, context);
      case STATIC_SPARSE_PARTITIONED_HYPERGRAPH:
        return internal::partition<LargeKHypergraphTypeTraits>(hypergraph, context);
      case DYNAMIC_PARTITIONED_HYPERGRAPH:
        return internal::partition<DynamicHypergraphTypeTraits>(hypergraph, context);
      case NULLPTR_PARTITION:
        return mt_kahypar_partitioned_hypergraph_t { nullptr, NULLPTR_PARTITION };
    }
    return mt_kahypar_partitioned_hypergraph_t { nullptr, NULLPTR_PARTITION };
  }

}  // namespace mt_kahypar

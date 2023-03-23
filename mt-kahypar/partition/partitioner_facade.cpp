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
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/csv_output.h"
#include "mt-kahypar/io/sql_plottools_serializer.h"
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

    // Partition Hypergraph
    PartitionedHypergraph partitioned_hg =
      Partitioner<TypeTraits>::partition(hg, context);

    return mt_kahypar_partitioned_hypergraph_t {
      reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(
        new PartitionedHypergraph(std::move(partitioned_hg))), PartitionedHypergraph::TYPE };
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

  void PartitionerFacade::printPartitioningResults(const mt_kahypar_partitioned_hypergraph_t phg,
                                                   const Context& context,
                                                   const std::chrono::duration<double>& elapsed_seconds) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      case STATIC_PARTITIONED_GRAPH:
        io::printPartitioningResults(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      case DYNAMIC_PARTITIONED_GRAPH:
        io::printPartitioningResults(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      case STATIC_PARTITIONED_HYPERGRAPH:
        io::printPartitioningResults(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      case STATIC_SPARSE_PARTITIONED_HYPERGRAPH:
        io::printPartitioningResults(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      case DYNAMIC_PARTITIONED_HYPERGRAPH:
        io::printPartitioningResults(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      case NULLPTR_PARTITION: break;
    }
  }

  std::string PartitionerFacade::serializeCSV(const mt_kahypar_partitioned_hypergraph_t phg,
                                              const Context& context,
                                              const std::chrono::duration<double>& elapsed_seconds) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      case STATIC_PARTITIONED_GRAPH:
        return io::csv::serialize(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      case DYNAMIC_PARTITIONED_GRAPH:
        return io::csv::serialize(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
      case STATIC_PARTITIONED_HYPERGRAPH:
        return io::csv::serialize(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
      case STATIC_SPARSE_PARTITIONED_HYPERGRAPH:
        return io::csv::serialize(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
      case DYNAMIC_PARTITIONED_HYPERGRAPH:
        return io::csv::serialize(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
      case NULLPTR_PARTITION: return "";
    }
    return "";
  }

  std::string PartitionerFacade::serializeResultLine(const mt_kahypar_partitioned_hypergraph_t phg,
                                                     const Context& context,
                                                     const std::chrono::duration<double>& elapsed_seconds) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      case STATIC_PARTITIONED_GRAPH:
        return io::serializer::serialize(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
      case DYNAMIC_PARTITIONED_GRAPH:
        return io::serializer::serialize(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
      case STATIC_PARTITIONED_HYPERGRAPH:
        return io::serializer::serialize(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
      case STATIC_SPARSE_PARTITIONED_HYPERGRAPH:
        return io::serializer::serialize(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
      case DYNAMIC_PARTITIONED_HYPERGRAPH:
        return io::serializer::serialize(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
      case NULLPTR_PARTITION: return "";
    }
    return "";
  }

  void PartitionerFacade::writePartitionFile(const mt_kahypar_partitioned_hypergraph_t phg,
                                             const std::string& filename) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      case STATIC_PARTITIONED_GRAPH:
        io::writePartitionFile(utils::cast_const<StaticPartitionedGraph>(phg), filename);
        break;
      case DYNAMIC_PARTITIONED_GRAPH:
        io::writePartitionFile(utils::cast_const<DynamicPartitionedGraph>(phg), filename);
        break;
      case STATIC_PARTITIONED_HYPERGRAPH:
        io::writePartitionFile(utils::cast_const<StaticPartitionedHypergraph>(phg), filename);
        break;
      case STATIC_SPARSE_PARTITIONED_HYPERGRAPH:
        io::writePartitionFile(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), filename);
        break;
      case DYNAMIC_PARTITIONED_HYPERGRAPH:
        io::writePartitionFile(utils::cast_const<DynamicPartitionedHypergraph>(phg), filename);
        break;
      case NULLPTR_PARTITION: break;
    }
  }

}  // namespace mt_kahypar

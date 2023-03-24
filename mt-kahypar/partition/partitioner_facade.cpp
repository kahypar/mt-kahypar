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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/csv_output.h"
#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/conversion.h"

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

} // namespace internal

  mt_kahypar_partitioned_hypergraph_t PartitionerFacade::partition(mt_kahypar_hypergraph_t hypergraph,
                                                                   Context& context) {
    const mt_kahypar_partition_type_t type = to_partition_c_type(
      context.partition.preset_type, context.partition.instance_type);
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return internal::partition<StaticGraphTypeTraits>(hypergraph, context);
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return internal::partition<DynamicGraphTypeTraits>(hypergraph, context);
      #endif
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return internal::partition<StaticHypergraphTypeTraits>(hypergraph, context);
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return internal::partition<LargeKHypergraphTypeTraits>(hypergraph, context);
      #endif
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return internal::partition<DynamicHypergraphTypeTraits>(hypergraph, context);
      #endif
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
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      #endif
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      #endif
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      #endif
      case NULLPTR_PARTITION: break;
    }
  }

  std::string PartitionerFacade::serializeCSV(const mt_kahypar_partitioned_hypergraph_t phg,
                                              const Context& context,
                                              const std::chrono::duration<double>& elapsed_seconds) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
      #endif
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return io::csv::serialize(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      case NULLPTR_PARTITION: return "";
    }
    return "";
  }

  std::string PartitionerFacade::serializeResultLine(const mt_kahypar_partitioned_hypergraph_t phg,
                                                     const Context& context,
                                                     const std::chrono::duration<double>& elapsed_seconds) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
      #endif
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      case NULLPTR_PARTITION: return "";
    }
    return "";
  }

  void PartitionerFacade::writePartitionFile(const mt_kahypar_partitioned_hypergraph_t phg,
                                             const std::string& filename) {
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<StaticPartitionedGraph>(phg), filename);
        break;
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<DynamicPartitionedGraph>(phg), filename);
        break;
      #endif
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<StaticPartitionedHypergraph>(phg), filename);
        break;
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        io::writePartitionFile(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), filename);
        break;
      #endif
      #ifdef KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<DynamicPartitionedHypergraph>(phg), filename);
        break;
      #endif
      case NULLPTR_PARTITION: break;
    }
  }

}  // namespace mt_kahypar

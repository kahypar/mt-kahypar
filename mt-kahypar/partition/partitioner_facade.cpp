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

#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/csv_output.h"
#include "mt-kahypar/io/sql_plottools_serializer.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {

namespace internal {

  template<typename TypeTraits>
  mt_kahypar_partitioned_hypergraph_t partition(mt_kahypar_hypergraph_t hypergraph,
                                                Context& context,
                                                TargetGraph* target_graph) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    Hypergraph& hg = utils::cast<Hypergraph>(hypergraph);

    // Partition Hypergraph
    PartitionedHypergraph partitioned_hg =
      Partitioner<TypeTraits>::partition(hg, context, target_graph);

    return mt_kahypar_partitioned_hypergraph_t {
      reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(
        new PartitionedHypergraph(std::move(partitioned_hg))), PartitionedHypergraph::TYPE };
  }

  template<typename TypeTraits>
  void improve(mt_kahypar_partitioned_hypergraph_t partitioned_hg,
               Context& context,
               TargetGraph* target_graph) {
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(partitioned_hg);

    // Improve partition
    Partitioner<TypeTraits>::partitionVCycle(phg, context, target_graph);
  }

  template<typename TypeTraits>
  mt_kahypar_partitioned_hypergraph_t create_partitioned_hypergraph(mt_kahypar_hypergraph_t hypergraph,
                                                                    const PartitionID k,
                                                                    const PartitionID* partition) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    Hypergraph& hg = utils::cast<Hypergraph>(hypergraph);

    PartitionedHypergraph partitioned_hg(k, hg, parallel_tag_t { });
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID& hn) {
      if ( partition[hn] < 0 || partition[hn] >= k ) {
        throw InvalidInputException(
          "Try to partition hypergraph into " + STR(k) + " blocks, "
          "but input block for node " + STR(hn) + " is: " + STR(partition[hn]));
      }
      partitioned_hg.setOnlyNodePart(hn, partition[hn]);
    });
    partitioned_hg.initializePartition();

    return mt_kahypar_partitioned_hypergraph_t {
      reinterpret_cast<mt_kahypar_partitioned_hypergraph_s*>(
        new PartitionedHypergraph(std::move(partitioned_hg))), PartitionedHypergraph::TYPE };
  }

  void check_if_feature_is_enabled(const mt_kahypar_partition_type_t type) {
    unused(type);
    #ifndef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
    if ( type == MULTILEVEL_GRAPH_PARTITIONING || type == N_LEVEL_GRAPH_PARTITIONING ) {
      throw InvalidParameterException(
        "Graph partitioning features are deactivated. Add -DKAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES=ON "
        "to the cmake command and rebuild Mt-KaHyPar.");
    }
    #endif
    #ifndef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
    if ( type == N_LEVEL_HYPERGRAPH_PARTITIONING || type == N_LEVEL_GRAPH_PARTITIONING ) {
      throw InvalidParameterException(
        "Quality preset features are deactivated. Add -DKAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES=ON "
        "to the cmake command and rebuild Mt-KaHyPar.");
    }
    #endif
    #ifndef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
    if ( type == LARGE_K_PARTITIONING ) {
      throw InvalidParameterException(
        "Large-k partitioning features are deactivated. Add -DKAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES=ON "
        "to the cmake command and rebuild Mt-KaHyPar.");
    }
    #endif
  }

  template<typename ReturnT, typename Func>
  ReturnT with_type_traits(mt_kahypar_partition_type_t type, Func f, const ReturnT& defaultVal) {
    internal::check_if_feature_is_enabled(type);
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return f(StaticGraphTypeTraits{});
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return f(StaticHypergraphTypeTraits{});
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return f(LargeKHypergraphTypeTraits{});
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return f(DynamicGraphTypeTraits{});
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return f(DynamicHypergraphTypeTraits{});
      #endif
      default:
        return defaultVal;
    }
    return defaultVal;
  }

  template<typename ReturnT, typename Func>
  ReturnT with_type_traits(mt_kahypar_hypergraph_t hypergraph, const Context& context, Func f, const ReturnT& defaultVal) {
    InstanceType hg_type = get_instance_type(hypergraph);
    ASSERT(hg_type == InstanceType::UNDEFINED
           || context.partition.instance_type == InstanceType::UNDEFINED
           || context.partition.instance_type == get_instance_type(hypergraph));

    InstanceType instance_type = hg_type == InstanceType::UNDEFINED ? context.partition.instance_type : get_instance_type(hypergraph);
    const mt_kahypar_partition_type_t type = to_partition_c_type(context.partition.preset_type, instance_type);
    return with_type_traits(type, f, defaultVal);
  }

  template<typename ReturnT, typename Func>
  ReturnT with_type_traits(mt_kahypar_partitioned_hypergraph_t phg, Func f, const ReturnT& defaultVal) {
    const mt_kahypar_partition_type_t type = phg.type;
    return with_type_traits(type, f, defaultVal);
  }

  template<typename Func>
  void with_type_traits(mt_kahypar_partitioned_hypergraph_t phg, Func f) {
    with_type_traits(phg, [f](auto typeTraits) {
      f(typeTraits);
      return 0; // dummy return that is ignored
    }, 0);
  }
} // namespace internal


  mt_kahypar_partitioned_hypergraph_t PartitionerFacade::partition(mt_kahypar_hypergraph_t hypergraph,
                                                                   Context& context,
                                                                   TargetGraph* target_graph) {
<<<<<<< HEAD
    mt_kahypar_partitioned_hypergraph_t defaultVal { nullptr, NULLPTR_PARTITION };

    return internal::with_type_traits(hypergraph, context, [&](auto typeTraits) {
      return internal::partition<decltype(typeTraits)>(hypergraph, context, target_graph);
    }, defaultVal);
=======
    const mt_kahypar_partition_type_t type = to_partition_c_type(
      context.partition.preset_type, context.partition.instance_type);
    internal::check_if_feature_is_enabled(type);
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return internal::partition<StaticGraphTypeTraits>(hypergraph, context, target_graph);
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return internal::partition<StaticHypergraphTypeTraits>(hypergraph, context, target_graph);
      case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return internal::partition<CompressedHypergraphTypeTraits>(hypergraph, context, target_graph);
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return internal::partition<LargeKHypergraphTypeTraits>(hypergraph, context, target_graph);
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return internal::partition<DynamicGraphTypeTraits>(hypergraph, context, target_graph);
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return internal::partition<DynamicHypergraphTypeTraits>(hypergraph, context, target_graph);
      #endif
      default:
        return mt_kahypar_partitioned_hypergraph_t { nullptr, NULLPTR_PARTITION };
    }
    return mt_kahypar_partitioned_hypergraph_t { nullptr, NULLPTR_PARTITION };
>>>>>>> origin/compression
  }

  void PartitionerFacade::improve(mt_kahypar_partitioned_hypergraph_t phg,
                                  Context& context,
                                  TargetGraph* target_graph) {
<<<<<<< HEAD
    internal::with_type_traits(phg, [&](auto typeTraits) {
      internal::improve<decltype(typeTraits)>(phg, context, target_graph);
    });
  }

  mt_kahypar_partitioned_hypergraph_t PartitionerFacade::createPartitionedHypergraph(mt_kahypar_hypergraph_t hypergraph,
                                                                                     const Context& context,
                                                                                     const PartitionID k,
                                                                                     const PartitionID* partition) {
    mt_kahypar_partitioned_hypergraph_t defaultVal { nullptr, NULLPTR_PARTITION };

    return internal::with_type_traits(hypergraph, context, [&](auto typeTraits) {
      return internal::create_partitioned_hypergraph<decltype(typeTraits)>(hypergraph, k, partition);
    }, defaultVal);
  }

  mt_kahypar_partitioned_hypergraph_t PartitionerFacade::readPartitionFile(mt_kahypar_hypergraph_t hypergraph,
                                                                           const Context& context,
                                                                           const std::string& filename) {
    const HypernodeID num_nodes = internal::with_type_traits(hypergraph, context, [&](auto typeTraits) {
      using HG = typename decltype(typeTraits)::Hypergraph;
      return utils::cast_const<HG>(hypergraph).initialNumNodes();
    }, 0);

    std::vector<PartitionID> partition;
    io::readPartitionFile(filename, num_nodes, partition);
    return createPartitionedHypergraph(hypergraph, context, context.partition.k, partition.data());
=======
    const mt_kahypar_partition_type_t type = to_partition_c_type(
      context.partition.preset_type, context.partition.instance_type);
    internal::check_if_feature_is_enabled(type);
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        internal::improve<StaticGraphTypeTraits>(partitioned_hg, context, target_graph); break;
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        internal::improve<StaticHypergraphTypeTraits>(partitioned_hg, context, target_graph); break;
      case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING:
        internal::improve<CompressedHypergraphTypeTraits>(partitioned_hg, context, target_graph); break;
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        internal::improve<LargeKHypergraphTypeTraits>(partitioned_hg, context, target_graph); break;
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        internal::improve<DynamicGraphTypeTraits>(partitioned_hg, context, target_graph); break;
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        internal::improve<DynamicHypergraphTypeTraits>(partitioned_hg, context, target_graph); break;
      #endif
      default: break;
    }
>>>>>>> origin/compression
  }

  void PartitionerFacade::printPartitioningResults(const mt_kahypar_partitioned_hypergraph_t phg,
                                                   const Context& context,
                                                   const std::chrono::duration<double>& elapsed_seconds) {
<<<<<<< HEAD
    internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      io::printPartitioningResults(utils::cast_const<PHG>(phg), context, elapsed_seconds);
    });
=======
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<CompressedPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        io::printPartitioningResults(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
        break;
      #endif
      default: break;
    }
>>>>>>> origin/compression
  }

  std::string PartitionerFacade::serializeCSV(const mt_kahypar_partitioned_hypergraph_t phg,
                                              const Context& context,
                                              const std::chrono::duration<double>& elapsed_seconds) {
<<<<<<< HEAD
    return internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      return io::csv::serialize(utils::cast_const<PHG>(phg), context, elapsed_seconds);
    }, std::string{});
=======
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
        break;
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
      case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<CompressedPartitionedHypergraph>(phg), context, elapsed_seconds);
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return io::csv::serialize(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return io::csv::serialize(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      default: return "";
    }
    return "";
>>>>>>> origin/compression
  }

  std::string PartitionerFacade::serializeResultLine(const mt_kahypar_partitioned_hypergraph_t phg,
                                                     const Context& context,
                                                     const std::chrono::duration<double>& elapsed_seconds) {
<<<<<<< HEAD
    return internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      return io::serializer::serialize(utils::cast_const<PHG>(phg), context, elapsed_seconds);
    }, std::string{});
=======
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<StaticPartitionedGraph>(phg), context, elapsed_seconds);
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<StaticPartitionedHypergraph>(phg), context, elapsed_seconds);
      case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<CompressedPartitionedHypergraph>(phg), context, elapsed_seconds);
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<DynamicPartitionedGraph>(phg), context, elapsed_seconds);
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        return io::serializer::serialize(utils::cast_const<DynamicPartitionedHypergraph>(phg), context, elapsed_seconds);
      #endif
      default: return "";
    }
    return "";
>>>>>>> origin/compression
  }

  void PartitionerFacade::writePartitionFile(const mt_kahypar_partitioned_hypergraph_t phg,
                                             const std::string& filename) {
<<<<<<< HEAD
    internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      io::writePartitionFile(utils::cast_const<PHG>(phg), filename);
    });
=======
    const mt_kahypar_partition_type_t type = phg.type;
    switch ( type ) {
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case MULTILEVEL_GRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<StaticPartitionedGraph>(phg), filename);
        break;
      #endif
      case MULTILEVEL_HYPERGRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<StaticPartitionedHypergraph>(phg), filename);
        break;
      case COMPRESSED_MULTILEVEL_HYPERGRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<CompressedPartitionedHypergraph>(phg), filename);
        break;
      #ifdef KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES
      case LARGE_K_PARTITIONING:
        io::writePartitionFile(utils::cast_const<StaticSparsePartitionedHypergraph>(phg), filename);
        break;
      #endif
      #ifdef KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES
      #ifdef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
      case N_LEVEL_GRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<DynamicPartitionedGraph>(phg), filename);
        break;
      #endif
      case N_LEVEL_HYPERGRAPH_PARTITIONING:
        io::writePartitionFile(utils::cast_const<DynamicPartitionedHypergraph>(phg), filename);
        break;
      #endif
      default: break;
    }
>>>>>>> origin/compression
  }

}  // namespace mt_kahypar

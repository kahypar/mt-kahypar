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
  mt_kahypar_partitioned_hypergraph_t read_partition_file(mt_kahypar_hypergraph_t hypergraph,
                                                          const Context& context,
                                                          const std::string& filename) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    Hypergraph& hg = utils::cast<Hypergraph>(hypergraph);

    std::vector<PartitionID> partition;
    io::readPartitionFile(filename, hg.initialNumNodes(), partition);

    PartitionedHypergraph partitioned_hg(context.partition.k, hg, parallel_tag_t { });
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID& hn) {
      partitioned_hg.setOnlyNodePart(hn, partition[hn]);
    });
    partitioned_hg.initializePartition();

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
  ReturnT with_type_traits(const Context& context, Func f, const ReturnT& defaultVal) {
    const mt_kahypar_partition_type_t type = to_partition_c_type(
      context.partition.preset_type, context.partition.instance_type);
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
    mt_kahypar_partitioned_hypergraph_t defaultVal { nullptr, NULLPTR_PARTITION };

    return internal::with_type_traits(context, [&](auto typeTraits) {
      return internal::partition<decltype(typeTraits)>(hypergraph, context, target_graph);
    }, defaultVal);
  }


  mt_kahypar_partitioned_hypergraph_t PartitionerFacade::readPartitionFile(mt_kahypar_hypergraph_t hypergraph,
                                                                           const Context& context,
                                                                           const std::string& filename) {
    mt_kahypar_partitioned_hypergraph_t defaultVal { nullptr, NULLPTR_PARTITION };

    return internal::with_type_traits(context, [&](auto typeTraits) {
      return internal::read_partition_file<decltype(typeTraits)>(hypergraph, context, filename);
    }, defaultVal);
  }


  void PartitionerFacade::improve(mt_kahypar_partitioned_hypergraph_t phg,
                                  Context& context,
                                  TargetGraph* target_graph) {
    internal::with_type_traits(phg, [&](auto typeTraits) {
      internal::improve<decltype(typeTraits)>(phg, context, target_graph);
    });
  }

  void PartitionerFacade::printPartitioningResults(const mt_kahypar_partitioned_hypergraph_t phg,
                                                   const Context& context,
                                                   const std::chrono::duration<double>& elapsed_seconds) {
    internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      io::printPartitioningResults(utils::cast_const<PHG>(phg), context, elapsed_seconds);
    });
  }

  std::string PartitionerFacade::serializeCSV(const mt_kahypar_partitioned_hypergraph_t phg,
                                              const Context& context,
                                              const std::chrono::duration<double>& elapsed_seconds) {
    return internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      return io::csv::serialize(utils::cast_const<PHG>(phg), context, elapsed_seconds);
    }, std::string{});
  }

  std::string PartitionerFacade::serializeResultLine(const mt_kahypar_partitioned_hypergraph_t phg,
                                                     const Context& context,
                                                     const std::chrono::duration<double>& elapsed_seconds) {
    return internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      return io::serializer::serialize(utils::cast_const<PHG>(phg), context, elapsed_seconds);
    }, std::string{});
  }

  void PartitionerFacade::writePartitionFile(const mt_kahypar_partitioned_hypergraph_t phg,
                                             const std::string& filename) {
    internal::with_type_traits(phg, [&](auto typeTraits) {
      using PHG = typename decltype(typeTraits)::PartitionedHypergraph;
      io::writePartitionFile(utils::cast_const<PHG>(phg), filename);
    });
  }

}  // namespace mt_kahypar

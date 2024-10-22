/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/utilities.h"


// yes this is hacky
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
#include "growt/allocator/alignedallocator.hpp"
#include "growt/data-structures/hash_table_mods.hpp"
#include "growt/data-structures/table_config.hpp"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

using hasher_type    = utils_tm::hash_tm::murmur2_hash;
using allocator_type = growt::AlignedAllocator<>;
using ConcurrentHashTable = typename growt::table_config<
  uint64_t, float, hasher_type, allocator_type, hmod::sync>::table_type;
using HashTableHandle = typename ConcurrentHashTable::handle_type;

namespace mt_kahypar {

using EdgeMetadata = float;

template<typename TypeTraits>
class Level {

  using Hypergraph = typename TypeTraits::Hypergraph;

public:
  explicit Level(Hypergraph&& contracted_hypergraph,
                 parallel::scalable_vector<HypernodeID>&& communities,
                 parallel::scalable_vector<EdgeMetadata>&& edge_md,
                 double coarsening_time) :
    _contracted_hypergraph(std::move(contracted_hypergraph)),
    _communities(std::move(communities)),
    _edge_metadata(std::move(edge_md)),
    _coarsening_time(coarsening_time) { }

  Hypergraph& contractedHypergraph() {
    return _contracted_hypergraph;
  }

  const Hypergraph& contractedHypergraph() const {
    return _contracted_hypergraph;
  }

  const parallel::scalable_vector<EdgeMetadata>& edgeMetadata() const {
    return _edge_metadata;
  }

  // ! Maps a global vertex id of the representative hypergraph
  // ! to its global vertex id in the contracted hypergraph
  HypernodeID mapToContractedHypergraph(const HypernodeID hn) const {
    ASSERT(hn < _communities.size());
    return _communities[hn];
  }

  double coarseningTime() const {
    return _coarsening_time;
  }

  void freeInternalData() {
    tbb::parallel_invoke([&] {
      _contracted_hypergraph.freeInternalData();
    }, [&] {
      parallel::free(_communities);
    });
  }

private:
  // ! Contracted Hypergraph
  Hypergraph _contracted_hypergraph;
  // ! Defines the communities that are contracted
  // ! in the coarse hypergraph
  parallel::scalable_vector<HypernodeID> _communities;
  // ! Metadata for guided coarsening (frequencies / ML results)
  parallel::scalable_vector<EdgeMetadata> _edge_metadata;
  // ! Time to create the coarsened hypergraph
  // ! (includes coarsening + contraction time)
  double _coarsening_time;
};

template<typename TypeTraits>
class UncoarseningData {

  using Hypergraph = typename TypeTraits::Hypergraph;
  using HypergraphFactory = typename Hypergraph::Factory;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using ParallelHyperedge = typename Hypergraph::ParallelHyperedge;

public:
  explicit UncoarseningData(bool n_level, Hypergraph& hg, parallel::scalable_vector<EdgeMetadata>&& edge_md, const Context& context) :
    nlevel(n_level),
    _hg(hg),
    _edge_metadata(std::move(edge_md)),
    _context(context) {
      if (n_level) {
        compactified_hg = std::make_unique<Hypergraph>();
        compactified_phg = std::make_unique<PartitionedHypergraph>();
      } else {
        size_t estimated_number_of_levels = UL(1);
        if ( hg.initialNumNodes() > context.coarsening.contraction_limit ) {
          estimated_number_of_levels = std::ceil( std::log2(
              static_cast<double>(hg.initialNumNodes()) /
              static_cast<double>(context.coarsening.contraction_limit)) /
            std::log2(context.coarsening.maximum_shrink_factor) ) + UL(1);
        }
        hierarchy.reserve(estimated_number_of_levels);
      }
      is_phg_initialized = false;
      partitioned_hg = std::make_unique<PartitionedHypergraph>();
    }

  ~UncoarseningData() noexcept {
    tbb::parallel_for(UL(0), hierarchy.size(), [&](const size_t i) {
      (hierarchy)[i].freeInternalData();
    }, tbb::static_partitioner());
  }

  void setPartitionedHypergraph(PartitionedHypergraph&& phg) {
    ASSERT(!is_phg_initialized);
    partitioned_hg = std::make_unique<PartitionedHypergraph>(std::move(phg));
    is_phg_initialized = true;
  }

  void finalizeCoarsening() {
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);
    if (nlevel) {
      // Create compactified hypergraph containing only enabled vertices and hyperedges
      // with consecutive IDs => Less complexity in initial partitioning.
      timer.start_timer("compactify_hypergraph", "Compactify Hypergraph");
      auto compactification = HypergraphFactory::compactify(_hg);
      *compactified_hg = std::move(compactification.first);
      compactified_hn_mapping = std::move(compactification.second);
      *compactified_phg = PartitionedHypergraph(_context.partition.k, *compactified_hg, parallel_tag_t());
      timer.stop_timer("compactify_hypergraph");
    } else {
      timer.start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
      // Free memory of temporary contraction buffer and
      // release coarsening memory in memory pool
      if (!hierarchy.empty()) {
        hierarchy.back().contractedHypergraph().freeTmpContractionBuffer();
      } else {
        _hg.freeTmpContractionBuffer();
      }
      if (_context.type == ContextType::main) {
        parallel::MemoryPool::instance().release_mem_group("Coarsening");
      }

      // Construct partitioned hypergraph for initial partitioning
      if ( !is_phg_initialized ) {
        *partitioned_hg = PartitionedHypergraph(_context.partition.k, _hg, parallel_tag_t());
      }
      if (!hierarchy.empty()) {
        partitioned_hg->setHypergraph(hierarchy.back().contractedHypergraph());
      }
      is_phg_initialized = true;
      timer.stop_timer("finalize_multilevel_hierarchy");
    }
    is_finalized = true;
  }

  void performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities, bool deterministic,
          const HighResClockTimepoint& round_start) {
    ASSERT(!is_finalized);
    Hypergraph& current_hg = hierarchy.empty() ? _hg : hierarchy.back().contractedHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities, deterministic);
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    vec<EdgeMetadata> contracted_md = accumulateMetadata(current_hg, contracted_hg, communities);
    hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), std::move(contracted_md), elapsed_time);
  }

  PartitionedHypergraph& coarsestPartitionedHypergraph() {
    if (nlevel) {
      return *compactified_phg;
    } else {
      return *partitioned_hg;
    }
  }

  const vec<EdgeMetadata>& coarsestEdgeMetadata() const {
    if (!hierarchy.empty()) {
      return hierarchy.back().edgeMetadata();
    } else {
      return _edge_metadata;
    }
  }

  // Multilevel Data
  vec<Level<TypeTraits>> hierarchy;

  // NLevel Data
  // ! Once coarsening terminates we generate a compactified hypergraph
  // ! containing only enabled vertices and hyperedges within a consecutive
  // ! ID range, which is then used for initial partitioning
  std::unique_ptr<Hypergraph> compactified_hg;
  // ! Mapping from vertex IDs of the original hypergraph to the IDs
  // ! in the compactified hypergraph
  vec<HypernodeID> compactified_hn_mapping;
  // ! Compactified partitioned hypergraph
  std::unique_ptr<PartitionedHypergraph> compactified_phg;
  // ! Contains timings how long a coarsening pass takes for each round
  vec<vec<ParallelHyperedge>> removed_hyperedges_batches;
  // ! Removed single-pin and parallel nets.
  // ! All hyperedges that are contained in one vector must be restored once
  // ! we completly processed a vector of batches.
  vec<double> round_coarsening_times;

  // Both
  bool is_phg_initialized;
  std::unique_ptr<PartitionedHypergraph> partitioned_hg;
  bool is_finalized = false;
  bool nlevel;

private:
  vec<EdgeMetadata> accumulateMetadata(const Hypergraph& current_hg, const Hypergraph& contracted_hg, const vec<HyperedgeID>& mapping) const {
    const vec<EdgeMetadata>& old_md = coarsestEdgeMetadata();
    vec<EdgeMetadata> new_md;
    if constexpr (Hypergraph::is_graph) {
      ConcurrentHashTable accumulator(2 * contracted_hg.initialNumEdges());
      tbb::parallel_invoke([&] {
        new_md.resize(contracted_hg.initialNumEdges());
      }, [&] {
        // accumulate the metadata
        current_hg.doParallelForAllEdges([&](HyperedgeID edge) {
          uint32_t source = mapping[current_hg.edgeSource(edge)];
          uint32_t target = mapping[current_hg.edgeTarget(edge)];
          if (source < target) {
            HashTableHandle handle = accumulator.get_handle();
            EdgeMetadata val = old_md[edge];
            uint64_t key = (static_cast<uint64_t>(source) << 32 | target);
            handle.insert_or_update(key, val, [=](float& old) { old += val; });
          }
        });
      });

      // write the new metadata into the vector
      contracted_hg.doParallelForAllEdges([&](HyperedgeID edge) {
        uint32_t source = contracted_hg.edgeSource(edge);
        uint32_t target = contracted_hg.edgeTarget(edge);
        if (source > target) {
          std::swap(source, target);
        }
        ALWAYS_ASSERT(source < target);
        uint64_t key = (static_cast<uint64_t>(source) << 32 | target);
        HashTableHandle handle = accumulator.get_handle();
        auto it = handle.find(key);
        ALWAYS_ASSERT(it != handle.end());
        auto result = *it;
        new_md[edge] = result.second;
      });
    }
    return new_md;
  }

  Hypergraph& _hg;
  vec<EdgeMetadata> _edge_metadata;
  const Context& _context;
};

typedef struct uncoarsening_data_s uncoarsening_data_t;

namespace uncoarsening {
  template<typename TypeTraits>
  uncoarsening_data_t* to_pointer(UncoarseningData<TypeTraits>& ip_data) {
    return reinterpret_cast<uncoarsening_data_t*>(&ip_data);
  }

  template<typename TypeTraits>
  UncoarseningData<TypeTraits>& to_reference(uncoarsening_data_t* ptr) {
    return *reinterpret_cast<UncoarseningData<TypeTraits>*>(ptr);
  }
}

}

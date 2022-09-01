/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/snodes_coarsening.h"

namespace mt_kahypar {
  using ds::SepNodesStack;
  using ds::SeparatedNodes;

class Level {

public:
  explicit Level(Hypergraph&& contracted_hypergraph,
                 parallel::scalable_vector<HypernodeID>&& communities,
                 double coarsening_time) :
    _contracted_hypergraph(std::move(contracted_hypergraph)),
    _communities(std::move(communities)),
    _coarsening_time(coarsening_time) { }

  Hypergraph& contractedHypergraph() {
    return _contracted_hypergraph;
  }

  const Hypergraph& contractedHypergraph() const {
    return _contracted_hypergraph;
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

  const vec<HypernodeID>& communities() const {
    return _communities;
  }

  vec<HypernodeID>& communities() {
    return _communities;
  }

private:
  // ! Contracted Hypergraph
  Hypergraph _contracted_hypergraph;
  // ! Defines the communities that are contracted
  // ! in the coarse hypergraph
  parallel::scalable_vector<HypernodeID> _communities;
  // ! Time to create the coarsened hypergraph
  // ! (includes coarsening + contraction time)
  double _coarsening_time;
};

class UncoarseningData {
public:
  explicit UncoarseningData(bool n_level, Hypergraph& hg, const Context& context) :
    nlevel(n_level),
    _hg(hg),
    _context(context) {
      if (n_level) {
        compactified_hg = std::make_unique<Hypergraph>();
        compactified_phg = std::make_unique<PartitionedHypergraph>();
      } else {
        size_t estimated_number_of_levels = 1UL;
        if ( hg.initialNumNodes() > context.coarsening.contraction_limit ) {
          estimated_number_of_levels = std::ceil( std::log2(
              static_cast<double>(hg.initialNumNodes()) /
              static_cast<double>(context.coarsening.contraction_limit)) /
            std::log2(context.coarsening.maximum_shrink_factor) ) + 1UL;
        }
        hierarchy.reserve(estimated_number_of_levels);
      }
      partitioned_hg = std::make_unique<PartitionedHypergraph>();
    }

  ~UncoarseningData() noexcept {
    tbb::parallel_for(0UL, hierarchy.size(), [&](const size_t i) {
      (hierarchy)[i].freeInternalData();
    }, tbb::static_partitioner());
  }

  void finalizeCoarsening() {
    if (nlevel) {
      // Create compactified hypergraph containing only enabled vertices and hyperedges
      // with consecutive IDs => Less complexity in initial partitioning.
      utils::Timer::instance().start_timer("compactify_hypergraph", "Compactify Hypergraph");
      auto compactification = HypergraphFactory::compactify(_hg);
      *compactified_hg = std::move(compactification.first);
      compactified_hn_mapping = std::move(compactification.second);
      *compactified_phg = PartitionedHypergraph(_context.partition.k, *compactified_hg, parallel_tag_t());
      utils::Timer::instance().stop_timer("compactify_hypergraph");
    } else {
      utils::Timer::instance().start_timer("finalize_multilevel_hierarchy", "Finalize Multilevel Hierarchy");
      // Free memory of temporary contraction buffer and
      // release coarsening memory in memory pool
      if (!hierarchy.empty()) {
        hierarchy.back().contractedHypergraph().freeTmpContractionBuffer();
      } else {
        _hg.freeTmpContractionBuffer();
      }
      if (_context.type == kahypar::ContextType::main) {
        parallel::MemoryPool::instance().release_mem_group("Coarsening");
      }

      Hypergraph& coarsest_hg = hierarchy.empty() ? _hg : hierarchy.back().contractedHypergraph();
      const HypernodeID separatedTargetSize = calculateSeparatedNodesTargetSize(coarsest_hg, _hg);
      if (separatedTargetSize < coarsest_hg.separatedNodes().onliest().numNodes()) {
        star_partitioning::coarsen(coarsest_hg.separatedNodes(),
                                   coarsest_hg, _context, separatedTargetSize);
      }
      if (_context.initial_partitioning.reinsert_separated
          && _context.type != kahypar::ContextType::main
          && coarsest_hg.separatedNodes().coarsest().numNodes() > 0) {
        // TODO(maas): might be problematic if hierarchy is empty
        hierarchy.emplace_back(HypergraphFactory::reinsertSeparatedNodes(coarsest_hg, coarsest_hg.separatedNodes().coarsest()),
                                                                         parallel::scalable_vector<HypernodeID>(), 0);
      }

      // Construct partitioned hypergraph for initial partitioning
      *partitioned_hg = PartitionedHypergraph(_context.partition.k, _hg.initialNumNodes() + _hg.numSeparatedNodes(),
                                              _hg.initialNumEdges() + 2 * _hg.numSeparatedEdges(), parallel_tag_t());
      partitioned_hg->setHypergraph(hierarchy.empty() ? _hg : hierarchy.back().contractedHypergraph());

      utils::Timer::instance().stop_timer("finalize_multilevel_hierarchy");
    }
    is_finalized = true;
  }

  void initializeRefinement() {
    if (_context.initial_partitioning.reinsert_separated
        && _context.type != kahypar::ContextType::main
        && !hierarchy.empty() && hierarchy.back().contractedHypergraph().numIncludedSeparated() > 0) {
      const Hypergraph& ip_graph = hierarchy.back().contractedHypergraph();
      ASSERT(&partitioned_hg->hypergraph() == &ip_graph);
      Hypergraph& hg = hierarchy.size() == 1 ? _hg : hierarchy[hierarchy.size() - 2].contractedHypergraph();
      ASSERT(ip_graph.firstIncludedSeparated() == hg.initialNumNodes());
      ds::Array<PartIdType> part_ids(ip_graph.initialNumNodes(), PartIdType(kInvalidPartition));
      partitioned_hg->extractPartIDs(part_ids);

      partitioned_hg->setHypergraph(hg);
      partitioned_hg->initializeSeparatedParts();
      const HypernodeID first_included_sep = ip_graph.firstIncludedSeparated();
      tbb::parallel_for(first_included_sep, ip_graph.initialNumNodes(), [&] (const HypernodeID& node) {
        partitioned_hg->separatedSetOnlyNodePart(node - first_included_sep, part_ids[node].load());
      });
      partitioned_hg->propagateSeparatedPartIDsToFinest();
      hierarchy.pop_back();
    }

    ASSERT([&] {
      for (HypernodeID node: partitioned_hg->nodes()) {
        if (partitioned_hg->partID(node) == kInvalidPartition) {
          return false;
        }
      }
      // if (partitioned_hg->hasSeparatedNodes()) {
      //   SeparatedNodes& separated_nodes = partitioned_hg->separatedNodes().finest();
      //   for (HypernodeID node = 0; node < separated_nodes.numNodes(); ++node) {
      //     if (partitioned_hg->separatedPartID(node) == kInvalidPartition) {
      //       return false;
      //     }
      //   }
      // }
      return true;
    }() );
  }

  void performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities,
          const HighResClockTimepoint& round_start) {
    ASSERT(!is_finalized);
    Hypergraph& current_hg = hierarchy.empty() ? _hg : hierarchy.back().contractedHypergraph();
    ASSERT(current_hg.initialNumNodes() == communities.size());
    Hypergraph contracted_hg = current_hg.contract(communities);
    const HighResClockTimepoint round_end = std::chrono::high_resolution_clock::now();
    const double elapsed_time = std::chrono::duration<double>(round_end - round_start).count();
    hierarchy.emplace_back(std::move(contracted_hg), std::move(communities), elapsed_time);
  }

  PartitionedHypergraph& coarsestPartitionedHypergraph() {
    if (nlevel) {
      return *compactified_phg;
    } else {
      return *partitioned_hg;
    }
  }

  HypernodeID calculateSeparatedNodesTargetSize(const Hypergraph& coarsened_hg, const Hypergraph& original_hg) {
    if (_context.coarsening.sep_nodes_coarsening_type == SNodesCoarseningSize::do_not_coarsen) {
      return coarsened_hg.numSeparatedNodes();
    }

    HypernodeID result = _context.coarsening.sep_nodes_coarsening_size_factor * coarsened_hg.initialNumNodes();
    // TODO: logarithmic does not work in initial partitioning
    if (_context.coarsening.sep_nodes_coarsening_type == SNodesCoarseningSize::logarithmic && original_hg.initialNumNodes() > 10000) {
      result *= std::log2(static_cast<double>(original_hg.initialNumNodes()) / 10000.0);
    }
    if (_context.type == kahypar::ContextType::main) {
      result *= _context.coarsening.sep_nodes_coarsening_relax_main_factor;
    }
    return std::max(result, 3 * _context.coarsening.contraction_limit_multiplier);
  }

  // Multilevel Data
  vec<Level> hierarchy;

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
  std::unique_ptr<PartitionedHypergraph> partitioned_hg;
  bool is_finalized = false;
  bool nlevel;

private:
  Hypergraph& _hg;
  const Context& _context;
};
}

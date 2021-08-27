/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
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

#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
namespace mt_kahypar {
class Level {

public:
  explicit Level(Hypergraph&& contracted_hypergraph,
                 parallel::scalable_vector<HypernodeID>&& communities,
                 double coarsening_time) :
    _representative_hypergraph(nullptr),
    _contracted_hypergraph(std::move(contracted_hypergraph)),
    _contracted_partitioned_hypergraph(),
    _communities(std::move(communities)),
    _coarsening_time(coarsening_time) { }

  void setRepresentativeHypergraph(PartitionedHypergraph* representative_hypergraph) {
    _representative_hypergraph = representative_hypergraph;
  }

  PartitionedHypergraph& representativeHypergraph() {
    ASSERT(_representative_hypergraph);
    return *_representative_hypergraph;
  }

  Hypergraph& contractedHypergraph() {
    return _contracted_hypergraph;
  }

  PartitionedHypergraph& contractedPartitionedHypergraph() {
    return _contracted_partitioned_hypergraph;
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
      _contracted_partitioned_hypergraph.freeInternalData();
    }, [&] {
      parallel::free(_communities);
    });
  }

private:
  // ! Hypergraph on the next finer level
  PartitionedHypergraph* _representative_hypergraph;
  // ! Contracted Hypergraph
  Hypergraph _contracted_hypergraph;
  // ! Partitioned Hypergraph
  PartitionedHypergraph _contracted_partitioned_hypergraph;
  // ! Defines the communities that are contracted
  // ! in the coarse hypergraph
  parallel::scalable_vector<HypernodeID> _communities;
  // ! Time to create the coarsened hypergraph
  // ! (includes coarsening + contraction time)
  double _coarsening_time;
};

class UncoarseningData {
public:
  explicit UncoarseningData(bool n_level, const Hypergraph& hg, const Context& context) :
    nlevel(n_level) {
      if (nlevel) {
        compactified_hg = std::make_shared<Hypergraph>();
        compactified_phg = std::make_shared<PartitionedHypergraph>();
        compactified_hn_mapping = std::make_shared<vec<HypernodeID>>();
        n_level_hierarchy = std::make_shared<VersionedBatchVector>();
        removed_hyperedges_batches = std::make_shared<vec<vec<ParallelHyperedge>>>();
        round_coarsening_times = std::make_shared<vec<double>>();
      } else {
        hierarchy = std::make_shared<vec<Level>>();
        size_t estimated_number_of_levels = 1UL;
        if ( hg.initialNumNodes() > context.coarsening.contraction_limit ) {
          estimated_number_of_levels = std::ceil( std::log2(
              static_cast<double>(hg.initialNumNodes()) /
              static_cast<double>(context.coarsening.contraction_limit)) /
            std::log2(context.coarsening.maximum_shrink_factor) ) + 1UL;
        }
        hierarchy->reserve(estimated_number_of_levels);
      }
      partitioned_hypergraph = std::make_shared<PartitionedHypergraph>();
    }

  // Multilevel Data
  std::shared_ptr<vec<Level>> hierarchy;

  // NLevel Data
  std::shared_ptr<Hypergraph> compactified_hg;
  std::shared_ptr<PartitionedHypergraph> compactified_phg;
  std::shared_ptr<vec<HypernodeID>> compactified_hn_mapping;
  std::shared_ptr<VersionedBatchVector> n_level_hierarchy;
  std::shared_ptr<vec<vec<ParallelHyperedge>>> removed_hyperedges_batches;
  std::shared_ptr<vec<double>> round_coarsening_times;

  // Both
  std::shared_ptr<PartitionedHypergraph> partitioned_hypergraph;
  bool is_finalized = false;
  bool nlevel;
};
}

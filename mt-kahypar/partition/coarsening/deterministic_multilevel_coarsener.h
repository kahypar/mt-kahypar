/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "multilevel_coarsener_base.h"
#include "i_coarsener.h"

#include "mt-kahypar/utils/reproducible_random.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/datastructures/buffered_vector.h"

#include <tbb/enumerable_thread_specific.h>

namespace mt_kahypar {
class DeterministicMultilevelCoarsener :  public ICoarsener,
                                          private MultilevelCoarsenerBase
{
public:
  DeterministicMultilevelCoarsener(Hypergraph& hypergraph,
                                   const Context& context,
                                   const bool top_level) :
    Base(hypergraph, context, top_level),
    propositions(hypergraph.initialNumNodes()),
    cluster_weight(hypergraph.initialNumNodes(), 0),
    opportunistic_cluster_weight(hypergraph.initialNumNodes(), 0),
    nodes_in_too_heavy_clusters(hypergraph.initialNumNodes()),
    default_rating_maps(hypergraph.initialNumNodes())
  {
  }

  ~DeterministicMultilevelCoarsener() {

  }

private:
  struct Proposition {
    HypernodeID node = kInvalidHypernode, cluster = kInvalidHypernode;
    HypernodeWeight weight = 0;
  };

  static constexpr bool debug = false;



  HypernodeID currentLevelContractionLimit() {
    const auto& hg = currentHypergraph();
    return std::max( _context.coarsening.contraction_limit,
               static_cast<HypernodeID>(
                    (hg.initialNumNodes() - hg.numRemovedHypernodes()) / _context.coarsening.maximum_shrink_factor) );
  }

  void coarsenImpl() override;

  void calculatePreferredTargetCluster(HypernodeID u, const vec<HypernodeID>& clusters);

  size_t approveVerticesInTooHeavyClusters(vec<HypernodeID>& clusters);

  PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm) override {
    return Base::doUncoarsen(label_propagation, fm);
  }

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
  }

  void setHierarchyImpl(std::shared_ptr<vec<Level>> hierarchy) override {
    _hierarchy = hierarchy;
  }

  void setPhgImpl(std::shared_ptr<PartitionedHypergraph> phg) override {
    _partitioned_hg = phg;
  }

  using Base = MultilevelCoarsenerBase;
  using Base::_context;

  utils::ParallelPermutation<HypernodeID> permutation;
  vec<HypernodeID> propositions;
  vec<HypernodeWeight> cluster_weight, opportunistic_cluster_weight;
  ds::BufferedVector<HypernodeID> nodes_in_too_heavy_clusters;
  tbb::enumerable_thread_specific<ds::SparseMap<HypernodeID, double>> default_rating_maps;
  tbb::enumerable_thread_specific<vec<HypernodeID>> ties;

};
}

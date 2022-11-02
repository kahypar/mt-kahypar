/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
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
                                   UncoarseningData& uncoarseningData) :
    Base(hypergraph, context, uncoarseningData),
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

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
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

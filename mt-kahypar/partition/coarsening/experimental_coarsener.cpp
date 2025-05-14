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

#include "experimental_coarsener.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

static constexpr bool debug = false;

template<typename TypeTraits>
bool ExperimentalCoarsener<TypeTraits>::coarseningPassImpl() {
  auto& timer = utils::Utilities::instance().getTimer(_context.utility_id);
  const auto pass_start_time = std::chrono::high_resolution_clock::now();
  timer.start_timer("coarsening_pass", "Clustering");

  // first, initialize the cluster ids
  const Hypergraph& hg = Base::currentHypergraph();
  DBG << V(_pass_nr)
      << V(hg.initialNumNodes())
      << V(hg.initialNumEdges())
      << V(hg.initialNumPins());

  size_t num_nodes = Base::currentNumNodes();
  const double num_nodes_before_pass = num_nodes;
  vec<HypernodeID> clusters(num_nodes, kInvalidHypernode);
  tbb::parallel_for(UL(0), num_nodes, [&](HypernodeID u) {
    // cluster_weight[u] = hg.nodeWeight(u);
    clusters[u] = u;
  });

  // START implementation of actual coarsening


  // END implementation of actual coarsening

  // This is how you can get randomization if you need it
  //
  // if ( _enable_randomization ) {
  //   utils::Randomize::instance().parallelShuffleVector( _current_vertices, UL(0), _current_vertices.size());
  // }

  // This assertion might be useful 
  //
  // HEAVY_COARSENING_ASSERT([&] {
  //   parallel::scalable_vector<HypernodeWeight> expected_weights(current_hg.initialNumNodes());
  //   // Verify that clustering is correct
  //   for ( const HypernodeID& hn : current_hg.nodes() ) {
  //     const HypernodeID u = hn;
  //     const HypernodeID root_u = cluster_ids[u];
  //     expected_weights[root_u] += current_hg.nodeWeight(hn);
  //   }

  //   // Verify that cluster weights are aggregated correct
  //   for ( const HypernodeID& hn : current_hg.nodes() ) {
  //     const HypernodeID u = hn;
  //     const HypernodeID root_u = cluster_ids[u];
  //     if ( root_u == u && expected_weights[u] != _cluster_weight[u] ) {
  //       LOG << "The expected weight of cluster" << u << "is" << expected_weights[u]
  //           << ", but currently it is" << _cluster_weight[u];
  //       return false;
  //     }
  //   }
  //   return true;
  // }(), "Clustering computed invalid cluster ids and weights");

  timer.stop_timer("coarsening_pass");
  ++_pass_nr;
  if (num_nodes_before_pass / num_nodes <= _context.coarsening.minimum_shrink_factor) {
    return false;
  }

  _timer.start_timer("contraction", "Contraction");
  // at this point, the coarsening is finished and we use the final cluster ids to perform the contraction
  _uncoarseningData.performMultilevelContraction(std::move(clusters), false /* deterministic */, pass_start_time);
  _timer.stop_timer("contraction");
  return true;
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(ExperimentalCoarsener)

}

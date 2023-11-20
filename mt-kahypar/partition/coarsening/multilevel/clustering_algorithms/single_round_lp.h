/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/parallel_pq.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/coarsening/multilevel/clustering_context.h"
#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/multilevel/num_nodes_tracker.h"


namespace mt_kahypar {

template<typename ScorePolicy, typename HeavyNodePenaltyPolicy, typename AcceptancePolicy>
class SingleRoundLP {
  using Rating = MultilevelVertexPairRater::Rating;
  using PQ = ds::MultiQueue<double, HypernodeID>;

 public:
  SingleRoundLP(const HypernodeID /*num_nodes*/, const Context& context):
    _context(context) { }

  template<bool has_fixed_vertices, typename Hypergraph, typename DegreeSimilarityPolicy>
  void performClustering(const Hypergraph& hg,
                         const parallel::scalable_vector<HypernodeID>& node_mapping,
                         const DegreeSimilarityPolicy& similarity_policy,
                         ClusteringContext<Hypergraph>& cc,
                         std::function<double (HypernodeID)> weight_ratio_for_node_fn = [](const HypernodeID) { return 1.0; },
                         int pass_nr = 0) {
    unused(pass_nr);
    unused(weight_ratio_for_node_fn);  // parameter only exists for compatibility with TwoHopClustering

    auto handle_node = [&](const HypernodeID hn) {
      // We perform rating if ...
      //  1.) The contraction limit of the current level is not reached
      //  2.) Vertex hn is not matched before
      if (hg.nodeIsEnabled(hn) && cc.shouldContinue() && cc.vertexIsUnmatched(hn)) {
        const Rating rating = cc.template rate<ScorePolicy, HeavyNodePenaltyPolicy,
                                               AcceptancePolicy, has_fixed_vertices>(hg, hn, similarity_policy);
        if (rating.target != kInvalidHypernode) {
          cc.template matchVertices<has_fixed_vertices>(hg, hn, rating.target);
        }
      }
    };

    if (_context.coarsening.prioritize_high_degree) {
      PQ parallel_pq(_context.shared_memory.num_threads);
      hg.doParallelForAllNodes([&](const HypernodeID hn) {
        double rating = 0;
        if (_context.coarsening.prioritize_with_edge_weight) {
          for (const HyperedgeID& he : hg.incidentEdges(hn)) {
            rating += hg.edgeWeight(he);
          }
        } else {
          rating = hg.nodeDegree(hn);
        }
        if (_context.coarsening.prioritize_with_node_weight) {
          rating /= hg.nodeWeight(hn);
        }
        parallel_pq.insert(rating, hn);
      });

      auto task = [&]{
        while (true) {
          HypernodeID hn;
          bool success = parallel_pq.tryPop(hn);
          if (!success) {
            return;
          }
          handle_node(hn);
        }
      };

      tbb::task_group tg;
      for (size_t i = 0; i < _context.shared_memory.num_threads; ++i) { tg.run(task); }
      tg.wait();
    } else {
      // We iterate in parallel over all vertices of the hypergraph and compute its contraction partner.
      tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID id) {
        ASSERT(id < node_mapping.size());
        handle_node(node_mapping[id]);
      });
    }
  }

 private:
  const Context& _context;
};

}  // namespace mt_kahypar

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
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_vertex_pair_rater.h"
#include "mt-kahypar/partition/coarsening/multilevel/num_nodes_tracker.h"


namespace mt_kahypar {

template<typename ScorePolicy, typename HeavyNodePenaltyPolicy, typename AcceptancePolicy>
class SingleRoundLP {
  using Rating = MultilevelVertexPairRater::Rating;

 public:
  SingleRoundLP(const HypernodeID /*num_nodes*/, const Context& context):
    _context(context) { }

  template<bool has_fixed_vertices, typename Hypergraph, typename DegreeSimilarityPolicy>
  void performClustering(const Hypergraph& hg,
                         const parallel::scalable_vector<HypernodeID>& node_mapping,
                         const HypernodeID hierarchy_contraction_limit,
                         vec<HypernodeID>& cluster_ids,
                         MultilevelVertexPairRater& rater,
                         ConcurrentClusteringData& clustering_data,
                         NumNodesTracker& num_nodes_tracker,
                         ds::FixedVertexSupport<Hypergraph>& fixed_vertices,
                         const DegreeSimilarityPolicy& similarity_policy) {
    // We iterate in parallel over all vertices of the hypergraph and compute its contraction partner.
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID id) {
      ASSERT(id < node_mapping.size());
      const HypernodeID hn = node_mapping[id];
      // We perform rating if ...
      //  1.) The contraction limit of the current level is not reached
      //  2.) Vertex hn is not matched before
      if (hg.nodeIsEnabled(hn)
          && num_nodes_tracker.currentNumNodes() > hierarchy_contraction_limit
          && clustering_data.vertexIsUnmatched(hn)) {
        const Rating rating = rater.template rate<ScorePolicy, HeavyNodePenaltyPolicy, AcceptancePolicy, has_fixed_vertices>(
                                     hg, hn, cluster_ids, clustering_data.clusterWeight(), fixed_vertices,
                                     similarity_policy, _context.coarsening.max_allowed_node_weight);
        if (rating.target != kInvalidHypernode) {
          bool success = clustering_data.template matchVertices<has_fixed_vertices>(
            hg, hn, rating.target, cluster_ids, rater, fixed_vertices);
          if (success) {
            // update the number of nodes in a way that minimizes synchronization overhead
            num_nodes_tracker.subtractNode(_context.shared_memory.original_num_threads, hierarchy_contraction_limit);
          }
        }
      }
    });
  }

 private:
  const Context& _context;
};

}  // namespace mt_kahypar

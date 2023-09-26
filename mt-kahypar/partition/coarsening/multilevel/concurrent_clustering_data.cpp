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

#include "mt-kahypar/partition/coarsening/multilevel/concurrent_clustering_data.h"

#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

namespace mt_kahypar {

template<typename Hypergraph>
void ConcurrentClusteringData::initializeCoarseningPass(Hypergraph& current_hg,
                                                        parallel::scalable_vector<HypernodeID>& cluster_ids) {
  cluster_ids.resize(current_hg.initialNumNodes());
  tbb::parallel_for(ID(0), current_hg.initialNumNodes(), [&](const HypernodeID hn) {
    // Reset clustering
    _matching_state[hn] = STATE(MatchingState::UNMATCHED);
    _matching_partner[hn] = hn;
    cluster_ids[hn] = hn;
    if ( current_hg.nodeIsEnabled(hn) ) {
      _cluster_weight[hn] = current_hg.nodeWeight(hn);
    }
  });
}

template<bool has_fixed_vertices, typename Hypergraph>
bool ConcurrentClusteringData::matchVertices(const Hypergraph& hypergraph,
                                             const HypernodeID u,
                                             const HypernodeID v,
                                             parallel::scalable_vector<HypernodeID>& cluster_ids,
                                             MultilevelVertexPairRater& rater,
                                             ds::FixedVertexSupport<Hypergraph>& fixed_vertices) {
  ASSERT(u < hypergraph.initialNumNodes());
  ASSERT(v < hypergraph.initialNumNodes());
  uint8_t unmatched = STATE(MatchingState::UNMATCHED);
  uint8_t match_in_progress = STATE(MatchingState::MATCHING_IN_PROGRESS);

  // Indicates that u wants to join the cluster of v.
  // Will be important later for conflict resolution.
  bool success = false;
  const HypernodeWeight weight_u = hypergraph.nodeWeight(u);
  HypernodeWeight weight_v = _cluster_weight[v];
  if ( weight_u + weight_v <= _context.coarsening.max_allowed_node_weight ) {
    if ( _matching_state[u].compare_exchange_strong(unmatched, match_in_progress) ) {
      _matching_partner[u] = v;
      // Current thread gets "ownership" for vertex u. Only threads with "ownership"
      // can change the cluster id of a vertex.

      uint8_t matching_state_v = _matching_state[v].load();
      if ( matching_state_v == STATE(MatchingState::MATCHED) ) {
        // Vertex v is already matched and will not change it cluster id any more.
        // In that case, it is safe to set the cluster id of u to the cluster id of v.
        const HypernodeID rep = cluster_ids[v];
        ASSERT(_matching_state[rep] == STATE(MatchingState::MATCHED));
        success = joinCluster<has_fixed_vertices>(hypergraph, u, rep, cluster_ids, fixed_vertices);
      } else if ( _matching_state[v].compare_exchange_strong(unmatched, match_in_progress) ) {
        // Current thread has the "ownership" for u and v and can change the cluster id
        // of both vertices thread-safe.
        success = joinCluster<has_fixed_vertices>(hypergraph, u, v, cluster_ids, fixed_vertices);
        _matching_state[v] = STATE(MatchingState::MATCHED);
      } else {
        // State of v must be either MATCHING_IN_PROGRESS or an other thread changed the state
        // in the meantime to MATCHED. We have to wait until the state of v changed to
        // MATCHED or resolve the conflict if u is matched within a cyclic matching dependency

        // Conflict Resolution
        while ( _matching_state[v] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {

          // Check if current vertex is in a cyclic matching dependency
          HypernodeID cur_u = u;
          HypernodeID smallest_node_id_in_cycle = cur_u;
          while ( _matching_partner[cur_u] != u && _matching_partner[cur_u] != cur_u ) {
            cur_u = _matching_partner[cur_u];
            smallest_node_id_in_cycle = std::min(smallest_node_id_in_cycle, cur_u);
          }

          // Resolve cyclic matching dependency
          // Vertex with smallest id starts to resolve conflict
          const bool is_in_cyclic_dependency = _matching_partner[cur_u] == u;
          if ( is_in_cyclic_dependency && u == smallest_node_id_in_cycle) {
            success = joinCluster<has_fixed_vertices>(hypergraph, u, v, cluster_ids, fixed_vertices);
            _matching_state[v] = STATE(MatchingState::MATCHED);
          }
        }

        // If u is still in state MATCHING_IN_PROGRESS its matching partner v
        // must be matched in the meantime with an other vertex. Therefore,
        // we try to match u with the representative v's cluster.
        if ( _matching_state[u] == STATE(MatchingState::MATCHING_IN_PROGRESS) ) {
          ASSERT( _matching_state[v] == STATE(MatchingState::MATCHED) );
          const HypernodeID rep = cluster_ids[v];
          success = joinCluster<has_fixed_vertices>(hypergraph, u, rep, cluster_ids, fixed_vertices);
        }
      }
      rater.markAsMatched(u);
      rater.markAsMatched(v);
      _matching_partner[u] = u;
      _matching_state[u] = STATE(MatchingState::MATCHED);
    }
  }
  return success;
}

// ! Only for testing
template<typename Hypergraph>
bool ConcurrentClusteringData::verifyClustering(const Hypergraph& current_hg,
                                                const parallel::scalable_vector<HypernodeID>& cluster_ids) const {
  parallel::scalable_vector<HypernodeWeight> expected_weights(current_hg.initialNumNodes());
  // Verify that clustering is correct
  for ( const HypernodeID& hn : current_hg.nodes() ) {
    const HypernodeID u = hn;
    const HypernodeID root_u = cluster_ids[u];
    if ( root_u != cluster_ids[root_u] ) {
      LOG << "Hypernode" << u << "is part of cluster" << root_u << ", but cluster"
          << root_u << "is also part of cluster" << cluster_ids[root_u];
      return false;
    }
    expected_weights[root_u] += current_hg.nodeWeight(hn);
  }

  // Verify that cluster weights are aggregated correct
  for ( const HypernodeID& hn : current_hg.nodes() ) {
    const HypernodeID u = hn;
    const HypernodeID root_u = cluster_ids[u];
    if ( root_u == u && expected_weights[u] != _cluster_weight[u] ) {
      LOG << "The expected weight of cluster" << u << "is" << expected_weights[u]
          << ", but currently it is" << _cluster_weight[u];
      return false;
    }
  }
  return true;
}

template<bool has_fixed_vertices, typename Hypergraph>
bool ConcurrentClusteringData::joinCluster(const Hypergraph& hypergraph,
                                           const HypernodeID u,
                                           const HypernodeID rep,
                                           vec<HypernodeID>& cluster_ids,
                                           ds::FixedVertexSupport<Hypergraph>& fixed_vertices) {
  ASSERT(rep == cluster_ids[rep]);
  bool success = false;
  const HypernodeWeight weight_of_u = hypergraph.nodeWeight(u);
  const HypernodeWeight weight_of_rep = _cluster_weight[rep];
  bool cluster_join_operation_allowed =
    weight_of_u + weight_of_rep <= _context.coarsening.max_allowed_node_weight;
  if constexpr ( has_fixed_vertices ) {
    if ( cluster_join_operation_allowed ) {
      cluster_join_operation_allowed = fixed_vertices.contract(rep, u);
    }
  }
  if ( cluster_join_operation_allowed ) {
    cluster_ids[u] = rep;
    _cluster_weight[rep] += weight_of_u;
    success = true;
  }
  _matching_partner[u] = u;
  _matching_state[u] = STATE(MatchingState::MATCHED);
  return success;
}

namespace {
  #define INITIALIZE_PASS(X) void ConcurrentClusteringData::initializeCoarseningPass(                         \
                                      X& current_hg, parallel::scalable_vector<HypernodeID>& cluster_ids)

  #define MATCH_VERTICES(X) bool ConcurrentClusteringData::matchVertices<false>(const X& hypergraph, const HypernodeID u, \
                                   const HypernodeID v, parallel::scalable_vector<HypernodeID>& cluster_ids,              \
                                   MultilevelVertexPairRater& rater, ds::FixedVertexSupport<X>& fixed_vertices)
  #define MATCH_VERTICES_FIXED(X) bool ConcurrentClusteringData::matchVertices<true>(const X& hypergraph,              \
                                         const HypernodeID u, const HypernodeID v,                                     \
                                         parallel::scalable_vector<HypernodeID>& cluster_ids,                          \
                                         MultilevelVertexPairRater& rater, ds::FixedVertexSupport<X>& fixed_vertices)

  #define VERIFY_CLUSTERING(X) bool ConcurrentClusteringData::verifyClustering(const X& current_hg,           \
                                        const parallel::scalable_vector<HypernodeID>& cluster_ids) const
}

INSTANTIATE_FUNC_WITH_HYPERGRAPHS(INITIALIZE_PASS)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(MATCH_VERTICES)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(MATCH_VERTICES_FIXED)
INSTANTIATE_FUNC_WITH_HYPERGRAPHS(VERIFY_CLUSTERING)

}  // namespace mt_kahypar

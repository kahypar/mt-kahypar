/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/process_mapping/greedy_mapping.h"

#include <numeric>
#include <queue>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_bitset.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

namespace {

static constexpr bool debug = false;

struct PQElement {
  HyperedgeWeight rating;
  HypernodeID u;
};

bool operator<(const PQElement& lhs, const PQElement& rhs) {
  return lhs.rating < rhs.rating || (lhs.rating == rhs.rating && lhs.u < rhs.u);
}

bool operator>(const PQElement& lhs, const PQElement& rhs) {
  return lhs.rating > rhs.rating || (lhs.rating == rhs.rating && lhs.u > rhs.u);
}

using PQ = std::priority_queue<PQElement>;


HypernodeID get_node_with_minimum_weighted_degree(const ds::StaticGraph& graph) {
  vec<HypernodeID> min_nodes;
  HyperedgeWeight min_weighted_degree = std::numeric_limits<HypernodeWeight>::max();
  for ( const HypernodeID& hn : graph.nodes() ) {
    HyperedgeWeight weighted_degree = 0;
    for ( const HyperedgeID he : graph.incidentEdges(hn) ) {
      weighted_degree += graph.edgeWeight(he);
    }
    if ( weighted_degree < min_weighted_degree ) {
      min_nodes.clear();
      min_nodes.push_back(hn);
      min_weighted_degree = weighted_degree;
    } else if ( weighted_degree == min_weighted_degree ) {
      min_nodes.push_back(hn);
    }
  }
  ASSERT(min_nodes.size() > 0);
  return min_nodes.size() == 1 ? min_nodes[0] :
    min_nodes[utils::Randomize::instance().getRandomInt(
      0, static_cast<int>(min_nodes.size() - 1), SCHED_GETCPU)];
}

template<typename CommunicationHypergraph>
void compute_greedy_mapping(CommunicationHypergraph& communication_hg,
                            const ProcessGraph& process_graph,
                            const Context&,
                            const HypernodeID seed_node) {
  // For each node u, the ratings store weight of all incident hyperedges
  // that connect u to partial assignment
  vec<HyperedgeWeight> rating(communication_hg.initialNumNodes(), 0);
  vec<bool> visited_hes(communication_hg.initialNumEdges(), false);
  vec<bool> up_to_date_ratings(communication_hg.initialNumNodes(), true);
  vec<HypernodeID> nodes_to_update;
  // Marks unassigned processors
  ds::Bitset unassigned_processors(process_graph.numBlocks());
  ds::StaticBitset unassigned_processors_view(
    unassigned_processors.numBlocks(), unassigned_processors.data());
  PQ pq;
  auto assign = [&](const HypernodeID u,
                    const PartitionID process) {
    ASSERT(process != kInvalidPartition && process < communication_hg.k());
    ASSERT(unassigned_processors.isSet(process));
    communication_hg.setNodePart(u, process);
    up_to_date_ratings[u] = false; // This marks u as assigned
    unassigned_processors.unset(process); // This marks the process as assigned

    DBG << "Assign node" << u << "to process" << process;

    // Update ratings
    nodes_to_update.clear();
    for ( const HyperedgeID& he : communication_hg.incidentEdges(u) ) {
      if ( !visited_hes[he] ) {
        const HyperedgeWeight edge_weight = communication_hg.edgeWeight(he);
        for ( const HypernodeID& pin : communication_hg.pins(he) ) {
          rating[pin] += edge_weight;
          if ( up_to_date_ratings[pin] ) {
            nodes_to_update.push_back(pin);
            up_to_date_ratings[pin] = false;
          }
        }
        visited_hes[he] = true;
      }
    }

    // Update PQ
    for ( const HypernodeID& hn : nodes_to_update ) {
      pq.push(PQElement { rating[hn], hn });
      up_to_date_ratings[hn] = true;
    }
  };

  communication_hg.resetPartition();
  // Initialize unassigned processors
  for ( PartitionID block = 0; block < process_graph.numBlocks(); ++block ) {
    unassigned_processors.set(block);
  }
  // Assign seed node to process with minimum weighted degree
  assign(seed_node, get_node_with_minimum_weighted_degree(process_graph.graph()));

  HyperedgeWeight actual_objective = 0;
  vec<PartitionID> tie_breaking;
  vec<HyperedgeWeight> tmp_ratings(communication_hg.initialNumNodes(), 0);
  while ( !pq.empty() ) {
    const PQElement best = pq.top();
    const HypernodeID u = best.u;
    pq.pop();

    if ( !up_to_date_ratings[u] ) {
      continue;
    }

    ASSERT(communication_hg.partID(u) == kInvalidPartition);
    // Assign node with the strongest connection to the partial assignment
    // to the process that minimizes the process mapping objective function.
    for ( const HyperedgeID& he : communication_hg.incidentEdges(u) ) {
      ds::Bitset& connectivity_set = communication_hg.deepCopyOfConnectivitySet(he);
      ds::StaticBitset connectivity_set_view(
        connectivity_set.numBlocks(), connectivity_set.data());
      const HyperedgeWeight edge_weight = communication_hg.edgeWeight(he);
      const HyperedgeWeight distance_before = process_graph.distance(connectivity_set_view);
      for ( const PartitionID process : unassigned_processors_view ) {
        const bool was_already_set = connectivity_set.isSet(process);
        connectivity_set.set(process);
        const HyperedgeWeight distance_after = process_graph.distance(connectivity_set_view);
        tmp_ratings[process] += (distance_after - distance_before) * edge_weight;
        if ( !was_already_set ) {
          connectivity_set.unset(process);
        }
      }
    }

    // Determine processor that would result in the least increase of the
    // process mapping objective function.
    HyperedgeWeight best_rating = std::numeric_limits<HyperedgeWeight>::max();
    for ( const PartitionID process : unassigned_processors_view ) {
      if ( tmp_ratings[process] < best_rating ) {
        tie_breaking.clear();
        tie_breaking.push_back(process);
        best_rating = tmp_ratings[process];
      } else if ( tmp_ratings[process] == best_rating ) {
        tie_breaking.push_back(process);
      }
      tmp_ratings[process] = 0;
    }

    // Assign node to processor that results in the least increase of the objective function
    ASSERT(tie_breaking.size() > 0);
    const PartitionID best_process = tie_breaking.size() == 1 ? tie_breaking[0] :
      tie_breaking[utils::Randomize::instance().getRandomInt(
        0, static_cast<int>(tie_breaking.size() - 1), SCHED_GETCPU)];
    actual_objective += best_rating;
    assign(u, best_process);
  }
  ASSERT(actual_objective == metrics::quality(communication_hg, Objective::process_mapping));
  ASSERT([&] {
    for ( const HypernodeID hn : communication_hg.nodes() ) {
      if ( communication_hg.partID(hn) == kInvalidHypernode ) {
        return false;
      }
    }
    return true;
  }(), "There are unassigned nodes");
  DBG << "Greedy mapping algorithm with seed node" << seed_node
      << "produced an mapping with solution quality" << actual_objective;
}

} // namespace

template<typename CommunicationHypergraph>
void GreedyMapping<CommunicationHypergraph>::mapToProcessGraph(CommunicationHypergraph& communication_hg,
                                                               const ProcessGraph& process_graph,
                                                               const Context& context) {
  ASSERT(communication_hg.initialNumNodes() == process_graph.graph().initialNumNodes());

  HyperedgeWeight best_objective = metrics::quality(communication_hg, Objective::process_mapping);
  vec<PartitionID> best_mapping(communication_hg.initialNumNodes(), 0);
  std::iota(best_mapping.begin(), best_mapping.end(), 0);
  for ( const HypernodeID& hn : communication_hg.nodes() ) {
    // Compute greedy mapping with the current node as seed node
    compute_greedy_mapping(communication_hg, process_graph, context, hn);

    // Check if new mapping is better than the currently best mapping
    const HyperedgeWeight objective = metrics::quality(communication_hg, Objective::process_mapping);
    if ( objective < best_objective ) {
      best_objective = objective;
      for ( const HypernodeID& u : communication_hg.nodes() ) {
        best_mapping[u] = communication_hg.partID(u);
      }
    }
  }

  // Apply best mapping
  communication_hg.resetPartition();
  for ( const HypernodeID& hn : communication_hg.nodes() ) {
    communication_hg.setOnlyNodePart(hn, best_mapping[hn]);
  }
  communication_hg.initializePartition();
}

INSTANTIATE_CLASS_WITH_PARTITIONED_HG(GreedyMapping)

}  // namespace kahypar
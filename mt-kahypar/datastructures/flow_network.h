#/*******************************************************************************
 * This file is part of KaHyPar.

 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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

#include <algorithm>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_array.h"
#include "external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h"
#include "external_tools/kahypar/kahypar/datastructure/sparse_set.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
namespace ds {

struct FlowEdge {
  NodeID source;
  NodeID target;
  Flow flow;
  Capacity capacity;
  size_t reverseEdge;

  void increaseFlow(const Flow delta_flow) {
    ASSERT(flow + delta_flow <= capacity, "Cannot increase flow above capacity!");
    flow += delta_flow;
  }
};

template <typename TypeTraits>
class FlowNetwork {
  using AdjacentList = std::vector<std::vector<FlowEdge> >;
  using ConstIncidenceIterator = std::vector<FlowEdge>::const_iterator;
  using IncidenceIterator = std::vector<FlowEdge>::iterator;
  using NodeIterator = std::pair<const NodeID*, const NodeID*>;
  using HypernodeIterator = std::pair<const HypernodeID*, const HypernodeID*>;
  using HyperGraph = typename TypeTraits::HyperGraph;

 public:
  static constexpr Flow kInfty = std::numeric_limits<Flow>::max() / 2;
  static constexpr NodeID kInvalidNode = std::numeric_limits<NodeID>::max() / 2;

  FlowNetwork(const HypernodeID initial_num_nodes, const HypernodeID initial_num_edges, const size_t size) :
    _initial_num_nodes(initial_num_nodes),
    _initial_num_edges(initial_num_edges),
    _initial_size(size),
    _num_nodes(0),
    _num_edges(0),
    _num_hyperedges(0),
    _num_undirected_edges(0),
    _total_weight_hyperedges(0),
    _nodes(size),
    _sources(size),
    _sinks(size),
    _hypernodes(initial_num_nodes),
    _removed_hypernodes(initial_num_nodes),
    _pins_block0(initial_num_edges, 0),
    _pins_block1(initial_num_edges, 0),
    _cur_block0(0),
    _cur_block1(1),
    _contains_graph_hyperedges(initial_num_nodes),
    _flow_graph(size, std::vector<FlowEdge>()),
    _visited(size),
    _he_visited(initial_num_edges) { }

  ~FlowNetwork() = default;

  FlowNetwork(const FlowNetwork&) = delete;
  FlowNetwork& operator= (const FlowNetwork&) = delete;

  FlowNetwork(FlowNetwork&&) = delete;
  FlowNetwork& operator= (FlowNetwork&&) = delete;

  // ################### Flow Network Stat ###################

  size_t numNodes() const {
    return _num_nodes;
  }

  size_t numEdges() const {
    return _num_edges;
  }

  size_t numUndirectedEdges() const {
    return _num_undirected_edges;
  }

  HyperedgeWeight totalWeightHyperedges() {
    return _total_weight_hyperedges;
  }

  size_t initialSize() const {
    return _initial_size;
  }

  // ################### Flow Network Construction ###################

  void buildFlowGraph(HyperGraph& hypergraph, const Context& context) {
    _visited.reset();
    for (const HypernodeID& hn : hypernodes()) {
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        if (!_visited[he]) {
          if (isHyperedgeOfSize(he, 1)) {
            addHyperedge(hypergraph, context, he);
            addPin(hypergraph, he, hn);
            _contains_graph_hyperedges.set(hn, true);
          } else if (hypergraph.edgeSize(he) == 2) {
            // Treat hyperedges of size 2 as graph edges
            // => remove hyperedge nodes in flow network
            addGraphEdge(hypergraph, he);
          } else {
            // Add a directed flow edge between the incomming
            // and outgoing hyperedge node of he with
            // capacity w(he)
            addHyperedge(hypergraph, context, he);
          }

          _visited.set(he, true);
        }
      }
    }

    for (const NodeID& node : nodes()) {
      if (interpreteHyperedge(node, false)) {
        const HyperedgeID he = mapToHyperedgeID(node);
        if (hypergraph.edgeSize(he) != 2 && !isHyperedgeOfSize(he, 1)) {
          // If degree of hypernode 'pin' is smaller than 3
          // add clique between all incident hyperedges of pin.
          // Otherwise connect pin with hyperedge in flow network.
          _he_visited.reset();
          for (const HypernodeID& pin : hypergraph.pins(he)) {
            if (_hypernodes.contains(pin)) {
              if (!_contains_graph_hyperedges[pin] &&
                  hypergraph.nodeDegree(pin) <= 3) {
                ASSERT(!containsNode(pin), "Pin " << pin << " of HE " << he
                                                  << " is already contained in flow problem!");
                addClique(hypergraph, he, pin);
              } else {
                addPin(hypergraph, he, pin);
                ASSERT(containsNode(pin), "Pin is not contained in flow problem!");
              }
            }
          }
          // Invariant: All outgoing edges of hyperedge 'he'
          //            are contained in the flow problem.
        }
      }
    }
  }

  HyperedgeWeight build(HyperGraph& hypergraph, const Context& context, const PartitionID block_0, const PartitionID block_1) {
    buildFlowGraph(hypergraph, context);

    ASSERT([&]() {
          for (const HypernodeID& hn : hypernodes()) {
            if (containsNode(hn) && _removed_hypernodes.contains(hn)) {
              LOG << "Flow network contains HN " << hn << " but it is marked as removed!";
              return false;
            } else if (!containsNode(hn) && !_removed_hypernodes.contains(hn)) {
              LOG << "Flow network not contains HN " << hn << " and it is not marked as removed!";
              return false;
            }
          }
          return true;
        } (), "Hypernodes not correctly added or removed from flow network!");

    // Add fixed vertices to source and sink
    /*
    for (const HypernodeID& hn : hypernodes()) {
      if (hypergraph.isFixedVertex(hn)) {
        if (containsNode(hn)) {
          if (hypergraph.partID(hn) == block_0) {
            addSource(hn);
          } else if (hypergraph.partID(hn) == block_1) {
            addSink(hn);
          }
        } else {
          for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
            if (hypergraph.partID(hn) == block_0) {
              ASSERT(containsNode(mapToIncommingHyperedgeID(he)), "Source is not contained in flow problem!");
              addSource(mapToIncommingHyperedgeID(he));
            } else if (hypergraph.partID(hn) == block_1) {
              ASSERT(containsNode(mapToOutgoingHyperedgeID(he)), "Sink is not contained in flow problem!");
              addSink(mapToOutgoingHyperedgeID(he));
            }
          }
        }
      }
    }*/

    const HyperedgeWeight cut = buildSourcesAndSinks(hypergraph, context, block_0, block_1);
    return cut;
  }

  void addHypernode(HyperGraph& hypergraph, const HypernodeID hn) {
    ASSERT(!containsHypernode(hn), "HN " << hn << " already contained in flow problem!");
    _hypernodes.add(hn);
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      if (hypergraph.edgeSize(he) == 2) {
        _contains_graph_hyperedges.set(hn, true);
      }
      if (hypergraph.partID(hn) == _cur_block0) {
        _pins_block0.update(he, 1);
      } else if (hypergraph.partID(hn) == _cur_block1) {
        _pins_block1.update(he, 1);
      }
    }
  }

  void removeHypernode(HyperGraph& hypergraph, const HypernodeID hn) {
    ASSERT(containsHypernode(hn));
    _hypernodes.remove(hn);
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      if (hypergraph.edgeSize(he) == 2) {
        _contains_graph_hyperedges.set(hn, false);
      }
      if (hypergraph.partID(hn) == _cur_block0) {
        _pins_block0.update(he, -1);
      } else if (hypergraph.partID(hn) == _cur_block1) {
        _pins_block1.update(he, -1);
      }
    }
  }

  void reset(const PartitionID block0 = 0, const PartitionID block1 = 1) {
    _num_nodes = 0;
    _num_edges = 0;
    _num_hyperedges = 0;
    _num_undirected_edges = 0;
    _total_weight_hyperedges = 0;
    _nodes.clear();
    _sources.clear();
    _sinks.clear();
    _hypernodes.clear();
    _removed_hypernodes.clear();
    _pins_block0.resetUsedEntries();
    _pins_block1.resetUsedEntries();
    _cur_block0 = block0;
    _cur_block1 = block1;
    _contains_graph_hyperedges.reset();
    _visited.reset();
  }

  bool isTrivialFlow() const {
    const size_t num_hyperedges_s_t = _sources.size() + _sinks.size();
    return num_hyperedges_s_t == 2 * _num_hyperedges;
  }

  // ################### Flow Network <-> Hypergraph Mapping ###################

  bool isHypernode(const NodeID node) const {
    return node < _initial_num_nodes;
  }

  bool containsHypernode(const HypernodeID hn) {
    return _hypernodes.contains(hn);
  }

  std::pair<const HypernodeID*, const HypernodeID*> hypernodes() {
    return std::make_pair(_hypernodes.begin(), _hypernodes.end());
  }

  bool interpreteHypernode(const NodeID node) const {
    return node < _initial_num_nodes && _nodes.contains(node);
  }

  NodeID mapToIncommingHyperedgeID(const HyperedgeID he) const {
    return _initial_num_nodes + he;
  }

  NodeID mapToOutgoingHyperedgeID(const HyperedgeID he) const {
    return _initial_num_nodes + _initial_num_edges + he;
  }

  HyperedgeID mapToHyperedgeID(const NodeID node) const {
    // ASSERT(node >= _hg.initialNumNodes(), "Node " << node << " is not a hyperedge node!");
    if (node >= _initial_num_nodes + _initial_num_edges) {
      return node - _initial_num_nodes - _initial_num_edges;
    } else {
      return node - _initial_num_nodes;
    }
    return node;
  }

  bool interpreteHyperedge(const NodeID node, bool outgoing = true) const {
    return (outgoing && node >= _initial_num_nodes + _initial_num_edges) ||
           (!outgoing && node >= _initial_num_nodes &&
            node < _initial_num_nodes + _initial_num_edges);
  }

  std::pair<const HypernodeID*, const HypernodeID*> removedHypernodes() const {
    return std::make_pair(_removed_hypernodes.begin(), _removed_hypernodes.end());
  }

  bool isRemovedHypernode(const HypernodeID hn) const {
    return _removed_hypernodes.contains(hn);
  }

  bool isHyperedgeOfSize(const HyperedgeID he, const size_t size) {
    return _pins_block0.get(he) + _pins_block1.get(he) == size;
  }

  size_t pinsNotInFlowProblem(const HyperGraph& hypergraph, const HyperedgeID he, const PartitionID block) {
    if (block == _cur_block0) {
      return hypergraph.pinCountInPart(he, _cur_block0) - _pins_block0.get(he);
    } else {
      return hypergraph.pinCountInPart(he, _cur_block1) - _pins_block1.get(he);
    }
  }

  // ################### Flow Network Nodes ###################

  NodeIterator nodes() {
    return std::make_pair(_nodes.begin(), _nodes.end());
  }

  void addNode(const NodeID node) {
    if (!containsNode(node)) {
      _nodes.add(node);
      _flow_graph[node].clear();
      _num_nodes++;
    }
  }

  bool containsNode(const NodeID node) {
    return _nodes.contains(node);
  }

  // ################### Flow Network Edges ###################

  Flow residualCapacity(const FlowEdge& e) {
    return e.capacity - e.flow;
  }

  void increaseFlow(FlowEdge& e, const Flow delta_flow) {
    e.increaseFlow(delta_flow);
    FlowEdge& revEdge = reverseEdge(e);
    revEdge.increaseFlow(-delta_flow);
  }

  FlowEdge & reverseEdge(const FlowEdge& e) {
    return _flow_graph[e.target][e.reverseEdge];
  }

  std::pair<IncidenceIterator, IncidenceIterator> incidentEdges(const NodeID u) {
    ASSERT(_nodes.contains(u), "Node " << u << " is not part of the flow graph!");
    return std::make_pair(_flow_graph[u].begin(), _flow_graph[u].end());
  }


  std::pair<ConstIncidenceIterator, ConstIncidenceIterator> incidentEdges(const NodeID u) const {
    ASSERT(_nodes.contains(u), "Node " << u << " is not part of the flow graph!");
    return std::make_pair(_flow_graph[u].cbegin(), _flow_graph[u].cend());
  }

  // ################### Source And Sink ###################

  void addSource(const NodeID u) {
    _sources.add(u);
  }

  void addSink(const NodeID u) {
    _sinks.add(u);
  }

  size_t numSources() {
    return _sources.size();
  }

  size_t numSinks() {
    return _sinks.size();
  }

  bool isSource(const NodeID node) {
    return _sources.contains(node);
  }

  bool isSink(const NodeID node) {
    return _sinks.contains(node);
  }

  NodeIterator sources() {
    return std::make_pair(_sources.begin(), _sources.end());
  }

  NodeIterator sinks() {
    return std::make_pair(_sinks.begin(), _sinks.end());
  }

 protected:
  HyperedgeWeight buildSourcesAndSinks(HyperGraph& hypergraph, const Context& context, const PartitionID block_0, const PartitionID block_1) {
    _visited.reset();
    HyperedgeWeight cut = 0;
    for (const HypernodeID& hn : hypernodes()) {
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        if (!_visited[he]) {
          const size_t pins_u_block0 = _pins_block0.get(he);
          const size_t pins_u_block1 = _pins_block1.get(he);
          const size_t pins_not_u_block0 = hypergraph.pinCountInPart(he, block_0) - pins_u_block0;
          const size_t pins_not_u_block1 = hypergraph.pinCountInPart(he, block_1) - pins_u_block1;
          const size_t connectivity = (pins_u_block0 + pins_not_u_block0 > 0) + (pins_u_block1 + pins_not_u_block1 > 0);

          if (context.partition.objective == kahypar::Objective::cut &&
              !isRemovableFromCut(hypergraph, he, block_0, block_1)) {
            // Case 1: Hyperedge he cannot be removed from cut
            //         of k-way partition.
            //         E.g., if he contains a block not equal to
            //         block_0 and block_1
            //         => add incoming hyperedge node as source
            //            and outgoing hyperedge node as sink
            ASSERT(containsNode(mapToIncommingHyperedgeID(he)), "Source is not contained in flow problem!");
            ASSERT(containsNode(mapToOutgoingHyperedgeID(he)), "Sink is not contained in flow problem!");
            addSource(mapToIncommingHyperedgeID(he));
            addSink(mapToOutgoingHyperedgeID(he));
          } else {
            // Hyperedge he is a cut hyperedge of the hypergraph.
            // => if he contains pins from block_0 not contained
            //   in the flow problem, we add the incoming hyperedge
            //   node as source
            // => if he contains pins from block_1 not contained
            //   in the flow problem, we add the outgoing hyperedge
            //   node as sink
            if (pins_not_u_block0 > 0) {
              ASSERT(containsNode(mapToIncommingHyperedgeID(he)), "Source is not contained in flow problem!");
              addSource(mapToIncommingHyperedgeID(he));
            }
            if (pins_not_u_block1 > 0) {
              ASSERT(containsNode(mapToOutgoingHyperedgeID(he)), "Sink is not contained in flow problem!");
              addSink(mapToOutgoingHyperedgeID(he));
            }
          }

          // Sum up the weight of all cut edges of the bipartition (block0,block1)
          if ((context.partition.objective == kahypar::Objective::km1 && connectivity > 1) ||
              (context.partition.objective == kahypar::Objective::cut && hypergraph.connectivity(he) > 1)) {
            cut += hypergraph.edgeWeight(he);
          }

          _visited.set(he, true);
        }
      }
    }

    return cut;
  }

  bool isRemovableFromCut(HyperGraph& hypergraph, const HyperedgeID he, const PartitionID block_0, const PartitionID block_1) {
    if (hypergraph.connectivity(he) > 2) {
      return false;
    }
    for (const PartitionID& part : hypergraph.connectivitySet(he)) {
      if (part != block_0 && part != block_1) {
        return false;
      }
    }
    return true;
  }

  void addEdge(const NodeID u, const NodeID v, const Capacity capacity, bool undirected = false) {
    addNode(u);
    addNode(v);
    FlowEdge e1;
    e1.source = u;
    e1.target = v;
    e1.flow = 0;
    e1.capacity = capacity;
    FlowEdge e2;
    e2.source = v;
    e2.target = u;
    e2.flow = 0;
    e2.capacity = (undirected ? capacity : 0);
    _flow_graph[u].push_back(e1);
    _flow_graph[v].push_back(e2);
    size_t e1_idx = _flow_graph[u].size() - 1;
    size_t e2_idx = _flow_graph[v].size() - 1;
    _flow_graph[u][e1_idx].reverseEdge = e2_idx;
    _flow_graph[v][e2_idx].reverseEdge = e1_idx;
    _num_edges += (undirected ? 2 : 1);
    _num_undirected_edges += (undirected ? 1 : 0);
    _total_weight_hyperedges += (capacity < kInfty ? capacity : 0);
  }

  /*
   * In following with denote with:
   *   - v -> a hypernode node
   *   - e' -> a incomming hyperedge node
   *   - e'' -> a outgoing hyperedge node
   */

  /*
   * Add a hyperedge to flow problem
   * Edge e' -> e''
   */
  void addHyperedge(HyperGraph& hypergraph, const Context& context, const HyperedgeID he, bool ignoreHESizeOne = false) {
    const NodeID u = mapToIncommingHyperedgeID(he);
    const NodeID v = mapToOutgoingHyperedgeID(he);
    if (!ignoreHESizeOne && isHyperedgeOfSize(he, 1)) {
      // Check if he is really a hyperedge of size 1 in flow problem
      ASSERT([&]() {
            size_t num_flow_hns = 0;
            for (const HypernodeID& pin : hypergraph.pins(he)) {
              if (containsHypernode(pin)) {
                num_flow_hns++;
              }
            }
            return (num_flow_hns == 1 ? true : false);
          } (), "Hyperedge " << he << " is not a size 1 hyperedge in flow problem!");

      if ((pinsNotInFlowProblem(hypergraph, he, _cur_block0) > 0 &&
           pinsNotInFlowProblem(hypergraph, he, _cur_block1) > 0) ||
          (context.partition.objective == kahypar::Objective::cut &&
           !isRemovableFromCut(hypergraph, he, _cur_block0, _cur_block1))) {
        addEdge(u, v, hypergraph.edgeWeight(he));
      } else {
        if (pinsNotInFlowProblem(hypergraph, he, _cur_block0) > 0) {
          addNode(u);
        } else if (pinsNotInFlowProblem(hypergraph, he, _cur_block1) > 0) {
          addNode(v);
        }
      }
    } else {
      addEdge(u, v, hypergraph.edgeWeight(he));
    }

    if (containsNode(u) || containsNode(v)) {
      _num_hyperedges++;
    }
  }

  void addGraphEdge(HyperGraph& hypergraph, const HyperedgeID he) {
    ASSERT(hypergraph.edgeSize(he) == 2, "Hyperedge " << he << " is not a graph edge!");

    NodeID u = kInvalidNode;
    NodeID v = kInvalidNode;
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      if (u == kInvalidNode) {
        u = pin;
      } else {
        v = pin;
      }
    }

    ASSERT(containsHypernode(u) || containsHypernode(v),
           "Hyperedge " << he << " is not in contained in the flow problem!");
    if (containsHypernode(v)) {
      std::swap(u, v);
    }
    addNode(u);

    if (containsHypernode(v)) {
      addEdge(u, v, hypergraph.edgeWeight(he), true);
    }
  }

  /*
   * Connecting a hyperedge with its pin
   * Edges e'' -> v and v -> e'
   */
  void addPin(HyperGraph& hypergraph, const HyperedgeID he, const HypernodeID pin) {
    const NodeID u = mapToIncommingHyperedgeID(he);
    const NodeID v = mapToOutgoingHyperedgeID(he);
    if (containsNode(u) && containsNode(v)) {
      addEdge(pin, u, kInfty);
      addEdge(v, pin, kInfty);
    } else {
      if (containsNode(u)) {
        ASSERT(_flow_graph[u].size() == 0, "Pin of size 1 hyperedge already added in flow graph!");
        addEdge(u, pin, hypergraph.edgeWeight(he));
      } else if (containsNode(v)) {
        ASSERT(_flow_graph[v].size() == 0, "Pin of size 1 hyperedge already added in flow graph!");
        addEdge(pin, v, hypergraph.edgeWeight(he));
      }
    }
    ASSERT(!(containsNode(u) || containsNode(v)) || containsNode(pin),
           "Pin " << pin << " should be added to flow graph!");
  }

  /*
   * Connecting the outgoing hyperedge node of he to all incomming
   * hyperedge nodes of incident hyperedges of hn.
   */
  void addClique(HyperGraph& hypergraph, const HyperedgeID he, const HypernodeID hn) {
    ASSERT(hypergraph.nodeDegree(hn) <= 3,
           "Hypernode " << hn << " should not be removed from flow problem!");
    _removed_hypernodes.add(hn);
    _he_visited.set(he, true);
    for (const HyperedgeID& target : hypergraph.incidentEdges(hn)) {
      if (!_he_visited[target]) {
        addEdge(mapToOutgoingHyperedgeID(he), mapToIncommingHyperedgeID(target), kInfty);
        _he_visited.set(target, true);
      }
    }
  }

  size_t _initial_num_nodes;
  size_t _initial_num_edges;
  size_t _initial_size;
  size_t _num_nodes;
  size_t _num_edges;
  size_t _num_hyperedges;
  size_t _num_undirected_edges;
  HyperedgeWeight _total_weight_hyperedges;

  kahypar::ds::SparseSet<NodeID> _nodes;
  kahypar::ds::SparseSet<NodeID> _sources;
  kahypar::ds::SparseSet<NodeID> _sinks;
  kahypar::ds::SparseSet<HypernodeID> _hypernodes;
  kahypar::ds::SparseSet<HypernodeID> _removed_hypernodes;
  kahypar::ds::FastResetArray<size_t> _pins_block0;
  kahypar::ds::FastResetArray<size_t> _pins_block1;
  PartitionID _cur_block0;
  PartitionID _cur_block1;
  kahypar::ds::FastResetFlagArray<> _contains_graph_hyperedges;
  AdjacentList _flow_graph;

  kahypar::ds::FastResetFlagArray<> _visited;
  kahypar::ds::FastResetFlagArray<> _he_visited;
};


}  // namespace ds
}  // namespace mt-kahypar

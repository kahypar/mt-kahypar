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
  size_t nextEdge; 

  void increaseFlow(const Flow delta_flow) {
    ASSERT(flow + delta_flow <= capacity, "Cannot increase flow above capacity!");
    flow += delta_flow;
  }
};

template<typename FlowNet>
class FlowEdgeIterator{

  private:
    FlowEdgeIterator(FlowNet & flow_network, size_t idx)
    :_flow_network(flow_network),
     _idx(idx) {
      if(_idx != FlowNet::kInvalidNode){
        _edge = &_flow_network.getEdge(_idx);
      }
     }

  public:
    static FlowEdgeIterator begin(FlowNet & flow_network, size_t idx){return FlowEdgeIterator(flow_network, idx);}
    static FlowEdgeIterator end(FlowNet & flow_network){return FlowEdgeIterator(flow_network, FlowNet::kInvalidNode);}
  
    FlowEdgeIterator& operator++() {
        _idx = _edge->nextEdge;
        if(_idx != FlowNet::kInvalidNode){
          _edge = &_flow_network.getEdge(_idx);
        }
        return *this;
    }

    bool operator==(FlowEdgeIterator e) const { return _idx == e._idx; }
    bool operator!=(FlowEdgeIterator e) const { return _idx != e._idx; }

    FlowEdge& operator*()  { return *_edge; }
    FlowEdge* operator->()  { return _edge; }

  private:
    FlowNet & _flow_network;
    FlowEdge* _edge;
    size_t _idx;
};

template <typename TypeTraits, typename FlowTypeTraits>
class FlowNetwork {
  using AdjacentList = parallel::scalable_vector<FlowEdge>;
  using ConstIncidenceIterator = parallel::scalable_vector<FlowEdge>::const_iterator;
  using IncidenceIterator = parallel::scalable_vector<FlowEdge>::iterator;
  using NodeIterator = std::pair<const NodeID*, const NodeID*>;
  using HypernodeIterator = std::pair<const HypernodeID*, const HypernodeID*>;
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using Derived = typename FlowTypeTraits::FlowNetwork;
  using FlowEdgeIter = FlowEdgeIterator<FlowNetwork<TypeTraits, FlowTypeTraits>>;

 public:
  static constexpr Flow kInfty = std::numeric_limits<Flow>::max() / 2;
  static constexpr NodeID kInvalidNode = std::numeric_limits<HypernodeID>::max() / 2;

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
    _flow_graph_size(0),
    _flow_graph(),
    _flow_graph_idx(0),
    _node_to_edge(size, kInvalidNode),
    _visited(size),
    _he_visited(initial_num_edges),
    _contains_aquired_node(initial_num_edges){
    }

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

  void buildFlowGraph(HyperGraph& hypergraph, const Context& context, const PartitionID block_0, Scheduler & scheduler) {
    _visited.reset();
    for (const HypernodeID& ogHn : hypernodes()) {
      const HypernodeID& hn = hypergraph.globalNodeID(ogHn);
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        if (!_visited[hypergraph.originalEdgeID(he)]) {
          if(true){
            //fix all pins of he's, that contain an already aquired pin (only OptScheduling)
            static_cast<Derived*>(this)->fixAquiredHyperEdge(hypergraph, block_0, scheduler, he);
          }
          if(!_contains_aquired_node[he]){
            if (isHyperedgeOfSize(hypergraph, he, 1)) {
              addHyperedge(hypergraph, context, he);
              addPin(hypergraph, he, hn);
              _contains_graph_hyperedges.set(ogHn, true);
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
          }

          _visited.set(hypergraph.originalEdgeID(he), true);
        }
      }
    }

    for (const NodeID& node : nodes()) {
      if (interpreteHyperedge(node, false)) {
        const HyperedgeID he_og = mapToHyperedgeID(node);
        const HyperedgeID he = hypergraph.globalEdgeID(he_og);
        if (hypergraph.edgeSize(he) != 2 && !isHyperedgeOfSize(hypergraph, he, 1) && !_contains_aquired_node[he]) {
          // If degree of hypernode 'pin' is smaller than 3
          // add clique between all incident hyperedges of pin.
          // Otherwise connect pin with hyperedge in flow network.
          _he_visited.reset();
          for (const HypernodeID& pin : hypergraph.pins(he)) {
            if (_hypernodes.contains(hypergraph.originalNodeID(pin))) {
              if (!_contains_graph_hyperedges[hypergraph.originalNodeID(pin)] &&
                  hypergraph.nodeDegree(pin) <= 3 && !containsNodeGlobalId(hypergraph, pin)) {
                ASSERT(!containsNodeGlobalId(hypergraph, pin), "Pin " << pin << " of HE " << he
                                                  << " is already contained in flow problem!");
                addClique(hypergraph, he, pin);
              } else {
                addPin(hypergraph, he, pin);
                ASSERT(containsNodeGlobalId(hypergraph, pin), "Pin is not contained in flow problem!");
              }
            }
          }
          // Invariant: All outgoing edges of hyperedge 'he'
          //            are contained in the flow problem.
        }
      }
    }
  }

  HyperedgeWeight build(HyperGraph& hypergraph, const Context& context, const PartitionID block_0,
    const PartitionID block_1, Scheduler & scheduler) {
    buildFlowGraph(hypergraph, context, block_0, scheduler);

    ASSERT([&]() {
          for (const HypernodeID& ogHn : hypernodes()) {
            const HypernodeID& hn = hypergraph.globalNodeID(ogHn);
            if (containsNodeGlobalId(hypergraph, hn) && _removed_hypernodes.contains(hypergraph.originalNodeID(hn))) {
              LOG << "Flow network contains HN " << hn << " but it is marked as removed!";
              return false;
            } else if (!containsNodeGlobalId(hypergraph, hn) && !_removed_hypernodes.contains(hypergraph.originalNodeID(hn))) {
              LOG << "Flow network not contains HN " << hn << " and it is not marked as removed!";
              return false;
            }
          }
          return true;
        } (), "Hypernodes not correctly added or removed from flow network!");


    const HyperedgeWeight cut = buildSourcesAndSinks(hypergraph, context, block_0, block_1);
    return cut;
  }

  std::vector<HypernodeWeight> get_aquired_part_weight(HyperGraph& hypergraph, size_t block_0,size_t block_1){
      std::vector<HypernodeWeight> part_weight(2, 0);
      for(auto hn: hypernodes()){
      size_t block = hypergraph.partID(hn);
      if(block == block_0){
          part_weight[0] += hypergraph.nodeWeight(hn);
      }else if(block == block_1){
          part_weight[1] += hypergraph.nodeWeight(hn);
      } else {
          ASSERT(false, "hypernode does not belong to block0 or block1");
      }
      }
      return part_weight;
  }

  void addHypernode(HyperGraph& hypergraph, const HypernodeID hn) {
    ASSERT(!containsHypernode(hypergraph, hn), "HN " << hn << " already contained in flow problem!");
    _hypernodes.add(hypergraph.originalNodeID(hn));
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      if (hypergraph.edgeSize(he) == 2) {
        _contains_graph_hyperedges.set(hypergraph.originalNodeID(hn), true);
      }
      if (hypergraph.partID(hn) == _cur_block0) {
        _pins_block0.update(hypergraph.originalEdgeID(he), 1);
      } else if (hypergraph.partID(hn) == _cur_block1) {
        _pins_block1.update(hypergraph.originalEdgeID(he), 1);
      }
    }
  }

  void removeHypernode(HyperGraph& hypergraph, const HypernodeID hn) {
    ASSERT(containsHypernode(hypergraph, hn));
    _hypernodes.remove(hypergraph.originalNodeID(hn));
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      if (hypergraph.edgeSize(he) == 2) {
        _contains_graph_hyperedges.set(hypergraph.originalNodeID(hn), false);
      }
      if (hypergraph.partID(hn) == _cur_block0) {
        _pins_block0.update(hypergraph.originalEdgeID(he), -1);
      } else if (hypergraph.partID(hn) == _cur_block1) {
        _pins_block1.update(hypergraph.originalEdgeID(he), -1);
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
    _contains_aquired_node.reset();
    _flow_graph_idx = 0;
    std::fill(_node_to_edge.begin(), _node_to_edge.end(), kInvalidNode);
  }

  void releaseHyperNodes(HyperGraph& hypergraph, PartitionID block_0, PartitionID block_1, Scheduler & scheduler){
    static_cast<Derived*>(this)->releaseHyperNodesImpl(hypergraph, block_0, block_1, scheduler);
  }

  bool isTrivialFlow() const {
    const size_t num_hyperedges_s_t = _sources.size() + _sinks.size();
    return num_hyperedges_s_t == 2 * _num_hyperedges;
  }

  // ################### Flow Network <-> Hypergraph Mapping ###################

  bool isHypernode(const NodeID node) const {
    return node < _initial_num_nodes;
  }

  bool containsHypernode(HyperGraph& hypergraph, const HypernodeID hn) {
    return _hypernodes.contains(hypergraph.originalNodeID(hn));
  }

  // _hypernodes contains the orginal HypernodesIDs, because the global HypernodeID is not continious in a multiple NUMA-Node Environment
  // to use the Hypernode IDs for operations on the Hypergraph, they have to be converted back:
  //
  // for (const HypernodeID& ogHn : flow_network.hypernodes()) {
  //    const HypernodeID& hn = hg.globalNodeID(ogHn);
  //
  std::pair<const HypernodeID*, const HypernodeID*> hypernodes() {
    return std::make_pair(_hypernodes.begin(), _hypernodes.end());
  }

  bool interpreteHypernode(const NodeID node) const {
    return node < _initial_num_nodes && _nodes.contains(node);
  }

  NodeID mapToIncommingHyperedgeID(HyperGraph& hypergraph, const HyperedgeID he) const {
    return _initial_num_nodes + hypergraph.originalEdgeID(he);
  }

  NodeID mapToOutgoingHyperedgeID(HyperGraph& hypergraph, const HyperedgeID he) const {
    return _initial_num_nodes + _initial_num_edges + hypergraph.originalEdgeID(he);
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


  bool isRemovedHypernode(HyperGraph& hypergraph, const HypernodeID hn) const {
    return _removed_hypernodes.contains(hypergraph.originalNodeID(hn));
  }

  bool isHyperedgeOfSize(HyperGraph& hypergraph, const HyperedgeID he, const size_t size) {
    return (_pins_block0.get(hypergraph.originalEdgeID(he)) + _pins_block1.get(hypergraph.originalEdgeID(he))) == size;
  }

  size_t pinsNotInFlowProblem(const HyperGraph& hypergraph, const HyperedgeID he, const PartitionID block) {
    if (block == _cur_block0) {
      return hypergraph.pinCountInPart(he, _cur_block0) - _pins_block0.get(hypergraph.originalEdgeID(he));
    } else {
      return hypergraph.pinCountInPart(he, _cur_block1) - _pins_block1.get(hypergraph.originalEdgeID(he));
    }
  }

  // ################### Flow Network Nodes ###################

  NodeIterator nodes() {
    return std::make_pair(_nodes.begin(), _nodes.end());
  }

  void addNodeId(const NodeID node) {
    if (!containsNodeId(node)) {
      _nodes.add(node);
      _num_nodes++;
    }
  }

  bool containsNodeGlobalId(HyperGraph& hypergraph, const HypernodeID node) {
    return _nodes.contains(hypergraph.originalNodeID(node));
  }

  bool containsNodeId(const NodeID id) {
    return _nodes.contains(id);
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
    return _flow_graph[e.reverseEdge];
  }

  std::pair<IncidenceIterator, IncidenceIterator> edges() {
    return std::make_pair(_flow_graph.begin(), _flow_graph.begin() + _flow_graph_idx);
  }

  std::pair<FlowEdgeIter, FlowEdgeIter> incidentEdges(NodeID node){
    return std::make_pair(FlowEdgeIter::begin(*this, _node_to_edge[node]), FlowEdgeIter::end(*this));
  }

  FlowEdge & getEdge(size_t idx){
    ASSERT(idx != kInvalidNode && idx < _flow_graph.size(), "Trying to get invalid Node!");
    return _flow_graph[idx];
  }

  /*
  size_t getFirstFlowEdge(NodeID node){
    return _node_to_edge[node];
  }

  FlowEdge & getEdge(size_t idx){
    ASSERT(idx != kInvalidNode && idx < _flow_graph.size(), "Trying to get invalid Node!");
    return _flow_graph[idx];
  }*/

  // ################### Source And Sink ###################

  //Calculate id accordingly before calling
  void addSourceWithId(const NodeID u) {
    _sources.add(u);
  }

  //Calculate id accordingly before calling
  void addSinkWithId(const NodeID u) {
    _sinks.add(u);
  }

  size_t numSources() {
    return _sources.size();
  }

  size_t numSinks() {
    return _sinks.size();
  }

  bool isGlobalIdSource(HyperGraph& hypergraph, const NodeID node) {
    return _sources.contains(hypergraph.originalNodeID(node));
  }

  bool isGlobalIdSink(HyperGraph& hypergraph, const NodeID node) {
    return _sinks.contains(hypergraph.originalNodeID(node));
  }

  bool isIdSource(const NodeID node) {
    return _sources.contains(node);
  }

  bool isIdSink(const NodeID node) {
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
    for (const HypernodeID& ogHn : hypernodes()) {
      const HypernodeID& hn = hypergraph.globalNodeID(ogHn);
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        if (!_visited[hypergraph.originalEdgeID(he)] && !_contains_aquired_node[he]) {
          const size_t pins_u_block0 = _pins_block0.get(hypergraph.originalEdgeID(he));
          const size_t pins_u_block1 = _pins_block1.get(hypergraph.originalEdgeID(he));
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
            ASSERT(containsNodeId(mapToIncommingHyperedgeID(hypergraph, he)), "Source is not contained in flow problem!");
            ASSERT(containsNodeId(mapToOutgoingHyperedgeID(hypergraph, he)), "Sink is not contained in flow problem!");
            addSourceWithId(mapToIncommingHyperedgeID(hypergraph, he));
            addSinkWithId(mapToOutgoingHyperedgeID(hypergraph, he));
          } else {
            // Hyperedge he is a cut hyperedge of the hypergraph.
            // => if he contains pins from block_0 not contained
            //   in the flow problem, we add the incoming hyperedge
            //   node as source
            // => if he contains pins from block_1 not contained
            //   in the flow problem, we add the outgoing hyperedge
            //   node as sink
            if (pins_not_u_block0 > 0) {
              ASSERT(containsNodeId(mapToIncommingHyperedgeID(hypergraph, he)), "Source is not contained in flow problem!"<< V(pins_not_u_block0)
              << V(pins_not_u_block1)<< V(pins_u_block0)<< V(pins_u_block1));
              addSourceWithId(mapToIncommingHyperedgeID(hypergraph, he));
            }
            if (pins_not_u_block1 > 0) {
              ASSERT(containsNodeId(mapToOutgoingHyperedgeID(hypergraph, he)), "Sink is not contained in flow problem!"<< V(pins_not_u_block0)
              << V(pins_not_u_block1)<< V(pins_u_block0)<< V(pins_u_block1));
              addSinkWithId(mapToOutgoingHyperedgeID(hypergraph, he));
            }
          }

          // Sum up the weight of all cut edges of the bipartition (block0,block1)
          if ((context.partition.objective == kahypar::Objective::km1 && connectivity > 1) ||
              (context.partition.objective == kahypar::Objective::cut && hypergraph.connectivity(he) > 1)) {
            cut += hypergraph.edgeWeight(he);
          }

          _visited.set(hypergraph.originalEdgeID(he), true);
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
    addNodeId(u);
    addNodeId(v);
    FlowEdge e1;
    e1.source = u;
    e1.target = v;
    e1.flow = 0;
    e1.capacity = capacity;
    e1.nextEdge = _node_to_edge[u];

    FlowEdge e2;
    e2.source = v;
    e2.target = u;
    e2.flow = 0;
    e2.capacity = (undirected ? capacity : 0);
    e2.nextEdge = _node_to_edge[v];

    e1.reverseEdge = _flow_graph_idx + 1;
    e2.reverseEdge = _flow_graph_idx;
 
    if(_flow_graph_idx + 1 < _flow_graph_size ){
      _flow_graph[_flow_graph_idx] = e1;
      _node_to_edge[u] = _flow_graph_idx++;
      _flow_graph[_flow_graph_idx] = e2;
      _node_to_edge[v] = _flow_graph_idx++;
    }else{
      _flow_graph.push_back(e1);
      _node_to_edge[u] = _flow_graph_idx++;
      _flow_graph.push_back(e2);
      _node_to_edge[v] = _flow_graph_idx++;
      _flow_graph_size += 2;
    }

    ASSERT(_flow_graph[_node_to_edge[u]].source == e1.source, "Inserttion not as expected!");
    ASSERT(_flow_graph[_node_to_edge[v]].source == e2.source, "Inserttion not as expected!");

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
    const NodeID u = mapToIncommingHyperedgeID(hypergraph, he);
    const NodeID v = mapToOutgoingHyperedgeID(hypergraph, he);
    if (!ignoreHESizeOne && isHyperedgeOfSize(hypergraph, he, 1)) {
      // Check if he is really a hyperedge of size 1 in flow problem
      ASSERT([&]() {
            size_t num_flow_hns = 0;
            for (const HypernodeID& pin : hypergraph.pins(he)) {
              if (containsHypernode(hypergraph, pin)) {
                num_flow_hns++;
              }
            }
            return (num_flow_hns == 1 ? true : false);
          } (), "Hyperedge " << he << " is not a size 1 hyperedge in flow problem!");

      if ((pinsNotInFlowProblem(hypergraph, he, _cur_block0) > 0 &&
           pinsNotInFlowProblem(hypergraph, he, _cur_block1) > 0) ||
          (context.partition.objective == kahypar::Objective::cut &&
           !isRemovableFromCut(hypergraph, he, _cur_block0, _cur_block1))||
           context.refinement.flow.algorithm == FlowAlgorithm::flow_opt) {
        addEdge(u, v, hypergraph.edgeWeight(he));
      } else {
        if (pinsNotInFlowProblem(hypergraph, he, _cur_block0) > 0) {
          addNodeId(u);
        } else if (pinsNotInFlowProblem(hypergraph, he, _cur_block1) > 0) {
          addNodeId(v);
        }
      }
    } else {
      addEdge(u, v, hypergraph.edgeWeight(he));
    }

    if (containsNodeId(u) || containsNodeId(v)) {
      _num_hyperedges++;
    }
  }

  void addGraphEdge(HyperGraph& hypergraph, const HyperedgeID he) {
    ASSERT(hypergraph.edgeSize(he) == 2, "Hyperedge " << he << " is not a graph edge!");

    HypernodeID u = kInvalidNode;
    HypernodeID v = kInvalidNode;
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      if (u == kInvalidNode) {
        u = pin;
      } else {
        v = pin;
      }
    }

    ASSERT(containsHypernode(hypergraph, u) || containsHypernode(hypergraph, v),
           "Hyperedge " << he << " is not in contained in the flow problem!");
    if (containsHypernode(hypergraph, v)) {
      std::swap(u, v);
    }
    addNodeId(hypergraph.originalNodeID(u));

    if (containsHypernode(hypergraph, v)) {
      addEdge(hypergraph.originalNodeID(u), hypergraph.originalNodeID(v), hypergraph.edgeWeight(he), true);
    }
  }

  /*
   * Connecting a hyperedge with its pin
   * Edges e'' -> v and v -> e'
   */
  void addPin(HyperGraph& hypergraph, const HyperedgeID he, const HypernodeID pin) {
    const NodeID u = mapToIncommingHyperedgeID(hypergraph, he);
    const NodeID v = mapToOutgoingHyperedgeID(hypergraph, he);
    if (containsNodeId(u) && containsNodeId(v)) {
      addEdge(hypergraph.originalNodeID(pin), u, kInfty);
      addEdge(v, hypergraph.originalNodeID(pin), kInfty);
    } else {
      if (containsNodeId(u)) {
        ASSERT(_node_to_edge[u] == kInvalidNode, "Pin of size 1 hyperedge already added in flow graph!");
        addEdge(u, hypergraph.originalNodeID(pin), hypergraph.edgeWeight(he));
      } else if (containsNodeId(v)) {
        ASSERT(_node_to_edge[v] == kInvalidNode, "Pin of size 1 hyperedge already added in flow graph!");
        addEdge(hypergraph.originalNodeID(pin), v, hypergraph.edgeWeight(he));
      }
    }
    ASSERT(!(containsNodeId(u) || containsNodeId(v)) || containsNodeGlobalId(hypergraph, pin),
           "Pin " << pin << " should be added to flow graph!" << V(containsNodeId(u)) << V(containsNodeId(v))
           << V(containsNodeGlobalId(hypergraph, pin)) << V(hypergraph.originalNodeID(pin)));
  }

  /*
   * Connecting the outgoing hyperedge node of he to all incomming
   * hyperedge nodes of incident hyperedges of hn.
   */
  void addClique(HyperGraph& hypergraph, const HyperedgeID he, const HypernodeID hn) {
    ASSERT(hypergraph.nodeDegree(hn) <= 3,
           "Hypernode " << hn << " should not be removed from flow problem!");
    _removed_hypernodes.add(hypergraph.originalNodeID(hn));
    _he_visited.set(hypergraph.originalEdgeID(he), true);
    for (const HyperedgeID& target : hypergraph.incidentEdges(hn)) {
      if (!_he_visited[hypergraph.originalEdgeID(target)]) {
        addEdge(mapToOutgoingHyperedgeID(hypergraph, he), mapToIncommingHyperedgeID(hypergraph, target), kInfty);
        _he_visited.set(hypergraph.originalEdgeID(target), true);
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
  size_t _flow_graph_size;
  AdjacentList _flow_graph;
  size_t _flow_graph_idx;
  
  parallel::scalable_vector<size_t> _node_to_edge;

  kahypar::ds::FastResetFlagArray<> _visited;
  kahypar::ds::FastResetFlagArray<> _he_visited;
  kahypar::ds::FastResetFlagArray<> _contains_aquired_node;

};

template <typename TypeTraits, typename FlowTypeTraits>
class MatchingFlowNetwork:public FlowNetwork<TypeTraits,FlowTypeTraits>  {
  using Base = FlowNetwork<TypeTraits,FlowTypeTraits>;
  using Base::Base;

  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using Scheduler = typename FlowTypeTraits::Scheduler;

  public:
    void releaseHyperNodesImpl(HyperGraph& hypergraph, PartitionID block_0, PartitionID block_1, Scheduler & scheduler){
      unused(hypergraph);
      unused(block_0);
      unused(block_1);
      unused(scheduler);
    }

    void fixAquiredHyperEdge(HyperGraph& hypergraph, const PartitionID block_0, Scheduler & scheduler, const HyperedgeID he){
      unused(hypergraph);
      unused(block_0);
      unused(scheduler);
      unused(he);
    }
};

template <typename TypeTraits, typename FlowTypeTraits>
class OptFlowNetwork:public FlowNetwork<TypeTraits,FlowTypeTraits>  {
  using Base = FlowNetwork<TypeTraits,FlowTypeTraits>;
  using Base::Base;

  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using Scheduler = typename FlowTypeTraits::Scheduler;

  public:
    void releaseHyperNodesImpl(HyperGraph& hypergraph, PartitionID block_0, PartitionID block_1, Scheduler & scheduler){
      std::vector<HypernodeWeight> aquired_part_weight = this->get_aquired_part_weight(hypergraph,block_0, block_1);
      scheduler.release_block_weight(block_0, block_1, aquired_part_weight[0]);
      scheduler.release_block_weight(block_1, block_0, aquired_part_weight[1]);
      for (const HypernodeID& hn : this->hypernodes()) {
        scheduler.releaseNode(hn);
      }  
    }

    void fixAquiredHyperEdge(HyperGraph& hypergraph, const PartitionID block_0, Scheduler & scheduler, const HyperedgeID he){
      for (const HypernodeID& pin : hypergraph.pins(he)) {
        if (!this->_hypernodes.contains(hypergraph.originalNodeID(pin)) && scheduler.isAquired(pin)) {
          this->_contains_aquired_node.set(he, true);
          fixNodes(hypergraph, he, block_0);
          break;
        }
      }
    }

  private:
    void fixNodes(HyperGraph& hypergraph, HyperedgeID he, const PartitionID block_0){
    const size_t pins_u_block0 = this->_pins_block0.get(hypergraph.originalEdgeID(he));
    const size_t pins_u_block1 = this->_pins_block1.get(hypergraph.originalEdgeID(he));

    const NodeID u = this->mapToIncommingHyperedgeID(hypergraph, he);
    const NodeID v = this->mapToOutgoingHyperedgeID(hypergraph, he);

    if(pins_u_block0 > 0){
      this->addNodeId(u);
      this->addSourceWithId(u);
    }
    if(pins_u_block1 > 0){
      this->addNodeId(v);
      this->addSinkWithId(v);
    }

    for (const HypernodeID& pin : hypergraph.pins(he)) {
      if (this->_hypernodes.contains(hypergraph.originalNodeID(pin))) {
        if(hypergraph.partID(pin) == block_0){
          this->addEdge(u, pin, this->kInfty);
        }else{
          this->addEdge(pin, v, this->kInfty);
        }
      }
    }
  }
};

}  // namespace ds
}  // namespace mt-kahypar

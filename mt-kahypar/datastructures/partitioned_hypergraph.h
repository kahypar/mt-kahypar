/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <atomic>
#include <type_traits>

#include "tbb/parallel_invoke.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

// Forward
template <typename Hypergraph,
          typename HypergraphFactory>
class NumaPartitionedHypergraph;

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory>
class PartitionedHypergraph {

  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");
  static_assert(!Hypergraph::is_numa_aware,  "Only non-numa-aware hypergraphs are allowed");

  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
  // ! Iterator to iterate over the set of communities contained in a hyperedge
  using CommunityIterator = typename Hypergraph::CommunityIterator;

  // ! Generic function that will be called if a hypernode v moves from a block from to a block to for
  // ! each incident net of the moved vertex v.
  // ! It will be called with the following arguments
  // !  1.) Hyperedge Weight
  // !  2.) Hyperedge Size
  // !  3.) Pin count in block from after move
  // !  4.) Pin count in block to after move
  // ! This function can be used to compute e.g. the delta in cut or km1 metric after a move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  using PinCountAtomic = parallel::IntegralAtomicWrapper<HypernodeID>;

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = false;
  static constexpr bool is_partitioned = true;

  explicit PartitionedHypergraph() :
    _k(0),
    _node(0),
    _hg(nullptr),
    connectivity_sets(0, 0) { }

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
    _k(k),
    _node(hypergraph.numaNode()),
    _hg(&hypergraph),
    part(hypergraph.initialNumNodes(), kInvalidPartition),
    pins_in_part(hypergraph.initialNumEdges() * k, PinCountAtomic(0)),
    connectivity_sets(hypergraph.initialNumEdges(), k) { }

  explicit PartitionedHypergraph(const PartitionID k,
                                 const TaskGroupID,
                                 Hypergraph& hypergraph) :
    _k(k),
    _node(hypergraph.numaNode()),
    _hg(&hypergraph),
    pins_in_part(),
    connectivity_sets(0, 0) {
    tbb::parallel_invoke([&] {
      part.resize(hypergraph.initialNumNodes(), kInvalidPartition);
    }, [&] {
      pins_in_part.assign(hypergraph.initialNumEdges() * k, PinCountAtomic(0));
    }, [&] {
      connectivity_sets = ConnectivitySets(hypergraph.initialNumEdges(), k);
    });
  }

  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other) :
    _k(other._k),
    _node(other._node),
    _hg(other._hg),
    part_weight(std::move(other.part_weight)),
    part(std::move(other.part)),
    pins_in_part(std::move(other.pins_in_part)),
    connectivity_sets(std::move(other.connectivity_sets)) { }

  PartitionedHypergraph & operator= (PartitionedHypergraph&& other) {
    _k = other._k;
    _node = other._node;
    _hg = other._hg;
    part_weight = std::move(other.part_weight),
    part = std::move(other.part),
    pins_in_part = std::move(other.pins_in_part);
    connectivity_sets = std::move(other.connectivity_sets);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

  // ! Returns the underlying hypergraph
  Hypergraph& hypergraph() {
    ASSERT(_hg);
    return *_hg;
  }

  void setHypergraph(Hypergraph& hypergraph) {
    _hg = &hypergraph;
  }

  // ! Number of NUMA hypergraphs
  size_t numNumaHypergraphs() const {
    return _hg->numNumaHypergraphs();
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _hg->initialNumNodes();
  }

  // ! Initial number of hypernodes on numa node
  HypernodeID initialNumNodes(const int node) const {
    return _hg->initialNumNodes(node);
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _hg->numRemovedHypernodes();
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _hg->initialNumEdges();
  }

  // ! Initial number of hyperedges on numa node
  HyperedgeID initialNumEdges(const int node) const {
    return _hg->initialNumEdges(node);
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _hg->initialNumPins();
  }

  // ! Initial number of pins on numa node
  HypernodeID initialNumPins(const int node) const {
    return _hg->initialNumPins(node);
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _hg->initialTotalVertexDegree();
  }

  // ! Initial sum of the degree of all vertices on numa node
  HypernodeID initialTotalVertexDegree(const int node) const {
    return _hg->initialTotalVertexDegree(node);
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _hg->totalWeight();
  }

  // ! Number of blocks this hypergraph is partitioned into
  PartitionID k() const {
    return _k;
  }

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const TaskGroupID task_group_id, const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllNodes(task_group_id, f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const TaskGroupID task_group_id, const F& f) const {
    _hg->doParallelForAllNodes(task_group_id, f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID task_group_id, const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllEdges(task_group_id, f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID task_group_id, const F& f) const {
    _hg->doParallelForAllEdges(task_group_id, f);
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return _hg->nodes();
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HypernodeIterator> nodes(const int node) const {
    return _hg->nodes(node);
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    return _hg->edges();
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HyperedgeIterator> edges(const int node) const {
    return _hg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    return _hg->incidentEdges(u);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    return _hg->pins(e);
  }

  // ! Returns a range to loop over the set of block ids contained in hyperedge e.
  IteratorRange<ConnectivitySets::Iterator> connectivitySet(const HyperedgeID e) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    return connectivity_sets.connectivitySet(local_id);
  }

  // ####################### Hypernode Information #######################

  // ! Returns for a vertex of the hypergraph its original vertex id
  // ! Can be used to map the global vertex ids to a consecutive range
  // ! of nodes between [0,|V|).
  HypernodeID originalNodeID(const HypernodeID u) const {
    return _hg->originalNodeID(u);
  }

  // ! Reverse operation of originalNodeID(u)
  HypernodeID globalNodeID(const HypernodeID u) const {
    return _hg->globalNodeID(u);
  }

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return _hg->nodeWeight(u);
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    const PartitionID block = partID(u);
    if ( block != kInvalidPartition ) {
      ASSERT(block < _k);
      const HypernodeWeight delta = weight - _hg->nodeWeight(u);
      part_weight[block] += delta;
    }
    _hg->setNodeWeight(u, weight);
  }


  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return _hg->nodeDegree(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return _hg->nodeIsEnabled(u);
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    _hg->enableHypernode(u);
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    _hg->disableHypernode(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Returns for a edge of the hypergraph its original edge id
  // ! Can be used to map the global edge ids to a consecutive range
  // ! of edges between [0,|E|).
  HypernodeID originalEdgeID(const HyperedgeID e) const {
    return _hg->originalEdgeID(e);
  }

  // ! Reverse operation of originalEdgeID(e)
  HypernodeID globalEdgeID(const HyperedgeID e) const {
    return _hg->globalEdgeID(e);
  }

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return _hg->edgeWeight(e);
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    _hg->setEdgeWeight(e, weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return _hg->edgeSize(e);
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return _hg->edgeIsEnabled(e);
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    _hg->enableHyperedge(e);
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    _hg->disableHyperedge(e);
  }

  // ####################### Uncontraction #######################

  void uncontract(const Memento&, parallel::scalable_vector<HyperedgeID>&) {
    ERROR("uncontract(memento,parallel_he) is not supported in partitioned hypergraph");
  }

  void uncontract(const std::vector<Memento>&,
                  parallel::scalable_vector<HyperedgeID>&,
                  const kahypar::ds::FastResetFlagArray<>&,
                  const bool) {
    ERROR("uncontract(...) is not supported in partitioned hypergraph");
  }

  void restoreDisabledHyperedgesThatBecomeNonParallel(
    const Memento&,
    parallel::scalable_vector<HyperedgeID>&,
    const kahypar::ds::FastResetFlagArray<>&) {
    ERROR("restoreDisabledHyperedgesThatBecomeNonParallel(...) is not supported"
          << "in partitioned hypergraph");
  }

  parallel::scalable_vector<HyperedgeID> findDisabledHyperedgesThatBecomeNonParallel(
    const Memento&,
    parallel::scalable_vector<HyperedgeID>&,
    const kahypar::ds::FastResetFlagArray<>&) {
    ERROR("findDisabledHyperedgesThatBecomeNonParallel(...) is not supported"
          << "in partitioned hypergraph");
    return parallel::scalable_vector<HyperedgeID>();
  }

  // ####################### Restore Hyperedges #######################

  // ! Restores an hyperedge of a certain size.
  void restoreEdge(const HyperedgeID he, const size_t size,
                   const HyperedgeID representative = kInvalidHyperedge) {
    _hg->restoreEdge(he, size, representative);
    for ( const HypernodeID& pin : pins(he) ) {
      incrementPinCountInPart(he, partID(pin));
    }
  }

  // ! Restores a single-pin hyperedge
  void restoreSinglePinHyperedge(const HyperedgeID he) {
    restoreEdge(he, 1);
  }

  void restoreParallelHyperedge(const HyperedgeID,
                                const Memento&,
                                parallel::scalable_vector<HyperedgeID>&,
                                const kahypar::ds::FastResetFlagArray<>* batch_hypernodes = nullptr) {
    unused(batch_hypernodes);
    ERROR("restoreParallelHyperedge(...) is not supported in partitioned hypergraph");
  }

  // ! Sets the block id of an unassigned vertex u and updates related partition information.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u, PartitionID id) {
    ASSERT(id != kInvalidPartition && id < _k && part[u] == kInvalidPartition);
    // TODO find a way to specify memory_order for add_and_fetch, since std::memory_order_relaxed suffices to detect a balance violation
    if (part_weight[id] += nodeWeight(u) <= max_part_weight[id]) {
      part[u] = id;
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        incrementPinCountInPart(he, id);
      }
      return true;
    } else {
      return false;
    }
  }

  // ! Sets the block id of an unassigned vertex u and updates related partition information.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u,
                   PartitionID id,
                   parallel::scalable_vector<PartitionedHypergraph>& hypergraphs) {
    ASSERT(id != kInvalidPartition && id < _k);
    if (part_weight[id] += nodeWeight(u) <= max_part_weight[id]) {
      part[u] = id;
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        common::hypergraph_of_edge(he, hypergraphs).incrementPinCountInPart(he, id);
      }
      return true;
    } else {
      return false;
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to', and updates related partition information.
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(Move& m, GlobalMoveTracker& moveTracker) {
    ASSERT(_hg->nodeIsEnabled(m.node), "Hypernode" << m.node << "is disabled");
    ASSERT(m.from != kInvalidPartition && m.from < _k);
    ASSERT(m.to != kInvalidPartition && m.to < _k);
    ASSERT(m.from != m.to);

    if (part_weight[m.to] += nodeWeight(m.node) <= max_part_weight[m.to]) {
      part[m.node] = m.to;
      part_weight[m.from] -= nodeWeight(m.node);

      const uint32_t move_id = moveTracker.insertMove(m);

      for ( const HyperedgeID& he : incidentEdges(m.node) ) {
        incrementPinCountInPart(he, m.to);
        decrementPinCountInPart(he, m.from);

        // update first move in
        std::atomic<MoveID>& fmi = first_move_in[he * _k + m.to];
        MoveID expected = fmi.load(std::memory_order_relaxed);
        while ((moveTracker.isIDStale(expected) || expected >= move_id) && !fmi.compare_exchange_weak(expected, move_id)) {  }

        // update last move out
        std::atomic<MoveID>& lmo = last_move_out[he * _k + m.from];
        expected = lmo.load(std::memory_order_relaxed);
        while (expected <= move_id && !lmo.compare_exchange_weak(expected, move_id)) { }
      }

      return true;
    } else {
      return false;
    }
  }

  uint32_t lastMoveOut(HyperedgeID he, PartitionID block) const {
    return last_move_out[he *_k + block].load(std::memory_order_relaxed);
  }

  uint32_t firstMoveIn(HyperedgeID he, PartitionID block) const {
    return first_move_in[he * _k + block].load(std::memory_order_relaxed);
  }

  HypernodeID remainingPinsFromBeginningOfMovePhase(HyperedgeID he, PartitionID block) const {
    return original_pins_minus_moved_out[he * _k + block].load(std::memory_order_relaxed);
  }

  void resetStoredMoveIDs() {
    for (auto& x : last_move_out) x.store(0, std::memory_order_relaxed);    // should be called very rarely
    for (auto& x : first_move_in) x.store(0, std::memory_order_relaxed);
  }

  void setRemainingOriginalPins() {
    // TODO try using memcpy
    for (size_t i = 0; i < pins_in_part.size(); ++i) {
      original_pins_minus_moved_out[i].store( pins_in_part[i].load(std::memory_order_relaxed), std::memory_order_relaxed );
    }
  }


  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u, PartitionID from, PartitionID to, const DeltaFunction& delta_func = NOOP_FUNC) {
    // TODO eliminate by rewriting LP refiner to not use DeltaFunction
    ASSERT(_hg->nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(from != kInvalidPartition && from < _k);
    ASSERT(to != kInvalidPartition && to < _k);
    ASSERT(from != to);

    if (part_weight[to] += nodeWeight(u) <= max_part_weight[to]) {
      part[u] = to;
      part_weight[from] -= nodeWeight(u);

      for ( const HyperedgeID& he : incidentEdges(u) ) {
        HypernodeID pin_count_after_to = incrementPinCountInPart(he, to);
        HypernodeID pin_count_after_from = decrementPinCountInPart(he, from);
        delta_func(he, edgeWeight(he), edgeSize(he), pin_count_after_from, pin_count_after_to);
      }
      return true;
    } else {
      return false;
    }
  }

  // ! Helper function to compute delta for cut-metric after changeNodePart
  static HyperedgeWeight cutDelta(const HyperedgeID,
                                  const HyperedgeWeight edge_weight,
                                  const HypernodeID edge_size,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after) {
    // TODO eliminate
    if ( edge_size > 1 ) {
      if (pin_count_in_to_part_after == edge_size) {
        return -edge_weight;
      } else if (pin_count_in_from_part_after == edge_size - 1 &&
                pin_count_in_to_part_after == 1) {
        return edge_weight;
      }
    }
    return 0;
  }

  // ! Helper function to compute delta for km1-metric after changeNodePart
  static HyperedgeWeight km1Delta(const HyperedgeID,
                                  const HyperedgeWeight edge_weight,
                                  const HypernodeID,
                                  const HypernodeID pin_count_in_from_part_after,
                                  const HypernodeID pin_count_in_to_part_after) {
    // TODO eliminate
    return (pin_count_in_to_part_after == 1 ? edge_weight : 0) +
           (pin_count_in_from_part_after == 0 ? -edge_weight : 0);
  }

  // ! Block which vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    return part[u];
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( connectivity(he) > 1 ) {
        return true;
      }
    }
    return false;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    return connectivity_sets.connectivity(local_id);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID id) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(id != kInvalidPartition && id < _k);
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    return pins_in_part[local_id * _k + id];
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID id) const {
    ASSERT(id != kInvalidPartition && id < _k);
    return part_weight[id];
  }

  HyperedgeWeight km1Gain(HypernodeID u, PartitionID p) const {
    return _affinity[u * _k + partID(u)].w1pins.load(std::memory_order_relaxed)
           - _affinity[u * _k + p].w0pins.load(std::memory_order_relaxed);
  }

  HyperedgeWeight affinity(HypernodeID u, PartitionID p) const {
    return _affinity[u * _k + p].w0pins.load(std::memory_order_relaxed);
  }

  // ! Reset partition (not thread-safe)    TODO does anyone use this?
  void resetPartition() {
    part.assign(part.size(), kInvalidPartition);
    for (auto& x : pins_in_part) x.store(0, std::memory_order_relaxed);  // TODO try using memcpy
    for (auto& x : part_weight) x.store(0, std::memory_order_relaxed);
    connectivity_sets.reset();
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);

    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    connectivity_sets.memoryConsumption(connectivity_set_node);

    // TODO finish this function when everything else is done

    parent->addChild("Part Info", sizeof(HypernodeWeight) * _k);
    //parent->addChild("Vertex Part Info", sizeof(VertexPartInfo) * _hg->initialNumNodes());
    parent->addChild("Pin Count In Part", sizeof(PinCountAtomic) * _k * _hg->initialNumEdges());
  }


  // TODO move this function somewhere else.
  // ####################### Extract Block #######################

  // ! Extracts a block of a partition as separate hypergraph.
  // ! It also returns a vertex-mapping from the original hypergraph to the sub-hypergraph.
  // ! If cut_net_splitting is activated, hyperedges that span more than one block (cut nets) are split, which is used for the connectivity metric.
  // ! Otherwise cut nets are discarded (cut metric).
  std::pair<Hypergraph, parallel::scalable_vector<HypernodeID> > extract(const TaskGroupID& task_group_id, PartitionID block, bool cut_net_splitting) {
    ASSERT(block != kInvalidPartition && block < _k);

    // Compactify vertex ids
    parallel::scalable_vector<HypernodeID> hn_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    parallel::scalable_vector<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
    HypernodeID num_hypernodes = 0;
    HypernodeID num_hyperedges = 0;
    tbb::parallel_invoke([&] {
      for ( const HypernodeID& hn : nodes() ) {
        if ( partID(hn) == block ) {
          hn_mapping[originalNodeID(hn)] = num_hypernodes++;
        }
      }
    }, [&] {
      for ( const HyperedgeID& he : edges() ) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          he_mapping[originalEdgeID(he)] = num_hyperedges++;
        }
      }
    });

    // Extract plain hypergraph data for corresponding block
    using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
    HyperedgeVector edge_vector;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weight.resize(num_hyperedges);
      doParallelForAllEdges(task_group_id, [&](const HyperedgeID he) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          const HyperedgeID original_id = originalEdgeID(he);
          hyperedge_weight[he_mapping[original_id]] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == block ) {
              edge_vector[he_mapping[original_id]].push_back(hn_mapping[originalNodeID(pin)]);
            }
          }
        }
      });
    }, [&] {
      hypernode_weight.resize(num_hypernodes);
      doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
        if ( partID(hn) == block ) {
          hypernode_weight[hn_mapping[originalNodeID(hn)]] = nodeWeight(hn);
        }
      });
    });

    // Construct hypergraph
    Hypergraph extracted_hypergraph = HypergraphFactory::construct(
      task_group_id, num_hypernodes, num_hyperedges,
      edge_vector, hyperedge_weight.data(), hypernode_weight.data());

    // Set community ids
    doParallelForAllNodes(task_group_id, [&](const HypernodeID& hn) {
      if ( partID(hn) == block ) {
        const HypernodeID extracted_hn =
          extracted_hypergraph.globalNodeID(hn_mapping[originalNodeID(hn)]);
        extracted_hypergraph.setCommunityID(extracted_hn, _hg->communityID(hn));
      }
    });
    extracted_hypergraph.initializeCommunities(task_group_id);
    if ( _hg->hasCommunityNodeMapping() ) {
      extracted_hypergraph.setCommunityNodeMapping(_hg->communityNodeMapping());
    }

    return std::make_pair(std::move(extracted_hypergraph), std::move(hn_mapping));
  }

 private:
  template <typename HyperGraph,
            typename HyperGraphFactory>
  friend class NumaPartitionedHypergraph;


  HypernodeID decrementPinCountInPart(const HyperedgeID e, const PartitionID p) {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeID local_hyperedge_id = common::get_local_position_of_edge(e);
    ASSERT(local_hyperedge_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    const HypernodeID pin_count_after = --pins_in_part[local_hyperedge_id * _k + p];
    if ( pin_count_after == 0 ) {
      // Connectivity of hyperedge decreased
      connectivity_sets.remove(local_hyperedge_id, p);
      for (HypernodeID u : pins(e)) {
        _affinity[u * _k + p].w0pins.fetch_sub(edgeWeight(e), std::memory_order_relaxed);
        _affinity[u * _k + p].w1pins.fetch_sub(edgeWeight(e), std::memory_order_relaxed);
      }
    } else if ( pin_count_after == 1 ) {
      for (HypernodeID u : pins(e)) {
        _affinity[u * _k + p].w1pins.fetch_add(edgeWeight(e), std::memory_order_relaxed);
      }
    }
    return pin_count_after;
  }

  HypernodeID incrementPinCountInPart(const HyperedgeID e, const PartitionID p) {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeID local_hyperedge_id = common::get_local_position_of_edge(e);
    ASSERT(local_hyperedge_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    const HypernodeID pin_count_after = ++pins_in_part[local_hyperedge_id * _k + p];
    if ( pin_count_after == 1 ) {
      // Connectivity of hyperedge increased
      connectivity_sets.add(local_hyperedge_id, p);
      for (HypernodeID u : pins(e)) {
        _affinity[u * _k + p].w0pins.fetch_sub(edgeWeight(e), std::memory_order_relaxed);
        _affinity[u * _k + p].w1pins.fetch_add(edgeWeight(e), std::memory_order_relaxed);
      }
    }
    return pin_count_after;
  }


  // ! Number of blocks
  PartitionID _k;

  // ! NUMA node of this partitioned hypergraph
  int _node;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg;

  // ! Weight and size information for all blocks.
  vec< std::atomic<HypernodeWeight> > part_weight;

  vec< HypernodeWeight > max_part_weight;

  // ! Current block IDs of the vertices
  vec< PartitionID > part;

  // ! For each hyperedge and each block, _pins_in_part stores the number of pins in that block
  vec< PinCountAtomic > pins_in_part;

  // ! For each hyperedge and block, the number of pins_in_part at the beginning of a move phase minus the number of moved out pins
  vec< PinCountAtomic > original_pins_minus_moved_out;    // TODO would like to get rid of this

  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block / the last move to remove a pin from that block
  vec< std::atomic<MoveID> > first_move_in, last_move_out;

  // TODO we probably don't need connectivity sets any more, except in the IP hypergraphs which don't need parallelism support
  // ! For each hyperedge, _connectivity_sets stores the set of blocks that the hyperedge spans
  ConnectivitySets connectivity_sets;


  struct AffinityInformation {
    std::atomic<HyperedgeWeight> w0pins, w1pins;
    AffinityInformation() : w0pins(0), w1pins(0) { }
  };
  // ! For each (vertex u, part i), the sum of edge weights for edges incident to u with zero/one pins in part i
  std::vector< AffinityInformation > _affinity;



};

} // namespace ds
} // namespace mt_kahypar
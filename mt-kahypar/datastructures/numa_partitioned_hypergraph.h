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

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          bool track_border_vertices = true>
class NumaPartitionedHypergraph {

  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");
  static_assert(Hypergraph::is_numa_aware,  "Only numa-aware hypergraphs are allowed");

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

  // Type Traits
  using UnderlyingHypergraph = typename Hypergraph::Hypergraph;
  using HardwareTopology = typename Hypergraph::HardwareTopology;
  using TBBNumaArena = typename Hypergraph::TBBNumaArena;
  using PartitionedHyperGraph = PartitionedHypergraph<UnderlyingHypergraph, track_border_vertices>;


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

  /*!
   * For each block \f$V_i\f$ of the \f$k\f$-way partition \f$\mathrm{\Pi} = \{V_1, \dots, V_k\}\f$,
   * a PartInfo object stores the number of hypernodes currently assigned to block \f$V_i\f$
   * as well as the sum of their weights.
   */
  class PartInfo {
   public:
    explicit PartInfo() :
      weight(0),
      size(0) { }

    bool operator== (const PartInfo& other) const {
      return weight == other.weight && size == other.size;
    }

    parallel::IntegralAtomicWrapper<HypernodeWeight> weight;
    parallel::IntegralAtomicWrapper<int64_t> size;
  };

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = true;
  static constexpr bool is_partitioned = true;

  explicit NumaPartitionedHypergraph() :
    _k(kInvalidPartition),
    _hg(nullptr),
    _part_info(),
    _hypergraphs() { }

  explicit NumaPartitionedHypergraph(const PartitionID k,
                                     Hypergraph& hypergraph) :
    _k(k),
    _hg(&hypergraph),
    _part_info(k),
    _hypergraphs() {
    for ( size_t node = 0; node < hypergraph.numNumaHypergraphs(); ++node) {
      _hypergraphs.emplace_back(k, hypergraph.numaHypergraph(node));
    }
  }

  explicit NumaPartitionedHypergraph(const PartitionID k,
                                     const TaskGroupID task_group_id,
                                     Hypergraph& hypergraph) :
    _k(k),
    _hg(&hypergraph),
    _part_info(k),
    _hypergraphs() {
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(task_group_id, [&](const int) {
          _hypergraphs.emplace_back();
        });

    // Construct NUMA partitioned hypergraphs in parallel
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node] = PartitionedHyperGraph(k, task_group_id, hypergraph.numaHypergraph(node));
        });
  }

  NumaPartitionedHypergraph(const NumaPartitionedHypergraph&) = delete;
  NumaPartitionedHypergraph & operator= (const NumaPartitionedHypergraph &) = delete;

  NumaPartitionedHypergraph(NumaPartitionedHypergraph&& other) :
    _k(other._k),
    _hg(other._hg),
    _part_info(std::move(other._part_info)),
    _hypergraphs(std::move(other._hypergraphs)) { }

  NumaPartitionedHypergraph & operator= (NumaPartitionedHypergraph&& other) {
    _k = other._k;
    _hg = other._hg;
    _part_info = std::move(other._part_info);
    _hypergraphs = std::move(other._hypergraphs);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

  // ! Returns the underlying hypergraph
  Hypergraph& hypergraph() {
    ASSERT(_hg);
    return *_hg;
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
    static_cast<const NumaPartitionedHypergraph&>(*this).doParallelForAllNodes(task_group_id, f);
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
    static_cast<const NumaPartitionedHypergraph&>(*this).doParallelForAllEdges(task_group_id, f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const TaskGroupID task_group_id, const F& f) const {
    _hg->doParallelForAllEdges(task_group_id, f);
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  ConcatenatedRange<IteratorRange<HypernodeIterator>> nodes() const {
    return _hg->nodes();
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph on a numa node
  IteratorRange<HypernodeIterator> nodes(const int node) const {
    return _hg->nodes(node);
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  ConcatenatedRange<IteratorRange<HyperedgeIterator>> edges() const {
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
    return hypergraph_of_edge(e).connectivitySet(e);
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
    _hg->setNodeWeight(u, weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return _hg->nodeDegree(u);
  }

  // ! Returns, if the corresponding vertex is high degree vertex
  bool isHighDegreeVertex(const HypernodeID u) const {
    return _hg->isHighDegreeVertex(u);
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
  }

  // ! Restores a single-pin hyperedge
  void restoreSinglePinHyperedge(const HyperedgeID he) {
    _hg->restoreSinglePinHyperedge(he);
  }

  void restoreParallelHyperedge(const HyperedgeID,
                                const Memento&,
                                parallel::scalable_vector<HyperedgeID>&,
                                const kahypar::ds::FastResetFlagArray<>* batch_hypernodes = nullptr) {
    unused(batch_hypernodes);
    ERROR("restoreParallelHyperedge(...) is not supported in partitioned hypergraph");
  }

  // ####################### Partition Information #######################

  /**
   * About moving a vertex:
   * The functions setNodePart(...) and changeNodePart(...) are implemented thread-safe
   * and lock-free. All involved data structures are based on atomics. We maintain
   * the following information when moving a node:
   *  1.) Block id of a vertex
   *  2.) Pin count in a block of a hyperedge
   *  3.) Connectivity of a hyperedge (number of blocks to which the pins of a hyperedge
   *      belongs to)
   *  4.) Connectivity set of a hyperedge (blocks to which the pins of a hyperedge belongs
   *      to)
   *  5.) Number of incident cut hyperedges of a vertex (number with hyperedges with
   *      connectivity greater than 1)
   *  6.) Block weights and sizes.
   *
   * In order, to change the block id of a vertex, we start by performing a CAS operation
   * on the block id of a vertex. This prevents that two threads are moving one vertex to
   * different blocks concurrently. If the operation succeeds we start updating all involved
   * data structures. Otherwise, we abort the operation and return false (indicating that
   * the move of the vertex failed). On succees, we update the local block weights of the
   * calling thread (see description below) and update 2.) - 5.) by iterating over each
   * hyperedge and increment or decrement the pin count of a block of a hyperedge. Both
   * operations return the pin count after the atomic update and based on that we can decide
   * if the hyperedge is still cut or not.
   *
   * It is guaranteed, that if all parallel vertex move operations are finished, the values
   * for 1.) - 5.) reflects the current state of the hypergraph. However, this is not the case
   * in a parallel setting. Since the update of all incident edges of the moved vertex is not
   * performed in a transactional/exclusive fashion, it can happen that threads that read one
   * of the values 2.) to 5.) have an intermediate view on those. Algorithms that
   * build on top of that have to consider that behavior and resolve or deal with that issue on
   * a higher layer.
   *
   * Another feature of changeNodePart is that someone can pass a generic delta function to compute
   * e.g. the delta in the cut- or km1-metric when moving vertices in parallel. E.g. one
   * can pass the following function to changeNodePart (as delta_func) to compute the delta
   * for the cut-metric of the move:
   *
   * HyperedgeWeight delta = 0;
   * auto cut_delta = [&](const HyperedgeWeight edge_weight,
   *                             const HypernodeID edge_size,
   *                             const HypernodeID pin_count_in_from_part_after,
   *                             const HypernodeID pin_count_in_to_part_after) {
   *   delta += mt_kahypar::Hypergraph::cutDelta(
   *     edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
   * };
   * hypergraph.changeNodePart(hn, from, to, cut_delta); // 'delta' will contain afterwards the delta
   *                                                     // for the cut metric
   *
   * Summing up all deltas of all parallel moves will be delta after the parallel execution phase.
   *
   */

  // ! Sets the block id of an unassigned vertex u.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u, const PartitionID id) {
    if ( hypergraph_of_vertex(u).setNodePart(u, id, _hypergraphs) ) {
      // Update block weights
      ++_part_info[id].size;
      _part_info[id].weight += nodeWeight(u);
      return true;
    } else {
      return false;
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    if ( hypergraph_of_vertex(u).changeNodePart(u, from, to, _hypergraphs, delta_func) ) {
      // Update block weights
      --_part_info[from].size;
      _part_info[from].weight -= nodeWeight(u);
      ++_part_info[to].size;
      _part_info[to].weight += nodeWeight(u);
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
    return (pin_count_in_to_part_after == 1 ? edge_weight : 0) +
           (pin_count_in_from_part_after == 0 ? -edge_weight : 0);
  }

  // ! Block which vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).partID(u);
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  TRUE_SPECIALIZATION(track_border_vertices, bool) isBorderNode(const HypernodeID u) const {
    return hypergraph_of_vertex(u).isBorderNode(u);
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  FALSE_SPECIALIZATION(track_border_vertices, bool) isBorderNode(const HypernodeID u) const {
    return hypergraph_of_vertex(u).isBorderNode(u, _hypergraphs);
  }

  // ! Number of incident cut hyperedges of vertex u
  TRUE_SPECIALIZATION(track_border_vertices, HyperedgeID) numIncidentCutHyperedges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).numIncidentCutHyperedges(u);
  }

  // ! Number of incident cut hyperedges of vertex u
  FALSE_SPECIALIZATION(track_border_vertices, HyperedgeID) numIncidentCutHyperedges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).numIncidentCutHyperedges(u, _hypergraphs);
  }

  // ! Initializes the number of cut hyperedges for each vertex
  // ! NOTE, this function have to be called after initial partitioning
  // ! and before local search.
  void initializeNumCutHyperedges() {
    for( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      _hypergraphs[node].initializeNumCutHyperedges(_hypergraphs);
    }
  }

  // ! Initializes the number of cut hyperedges for each vertex
  // ! NOTE, this function have to be called after initial partitioning
  // ! and before local search.
  void initializeNumCutHyperedges(const TaskGroupID task_group_id) {
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
          _hypergraphs[node].initializeNumCutHyperedges(task_group_id, _hypergraphs);
        });
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    return hypergraph_of_edge(e).connectivity(e);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID id) const {
    return hypergraph_of_edge(e).pinCountInPart(e, id);
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID id) const {
    ASSERT(id != kInvalidPartition && id < _k);
    return _part_info[id].weight;
  }

  // ! Number of vertices in a block
  size_t partSize(const PartitionID id) const {
    ASSERT(id != kInvalidPartition && id < _k);
    return _part_info[id].size;
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    for ( PartitionedHyperGraph& partitioned_hypergraph : _hypergraphs  ) {
      partitioned_hypergraph.resetPartition();
    }

    // Reset block weights and sizes
    for ( PartitionID block = 0; block < _k; ++block ) {
      _part_info[block].weight = 0;
      _part_info[block].size = 0;
    }
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Part Info", sizeof(PartInfo) * _k);
    for ( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      utils::MemoryTreeNode* numa_hypergraph_node = parent->addChild(
        "NUMA Partitioned Hypergraph " + std::to_string(node));
      _hypergraphs[node].memoryConsumption(numa_hypergraph_node);
    }
  }

 private:

  // ####################### Helper Functions #######################

  const PartitionedHyperGraph& hypergraph_of_vertex(const HypernodeID u) const {
    int node = common::get_numa_node_of_vertex(u);
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node];
  }

  PartitionedHyperGraph& hypergraph_of_vertex(const HypernodeID u) {
    return const_cast<PartitionedHyperGraph&>(static_cast<const NumaPartitionedHypergraph&>(*this).hypergraph_of_vertex(u));
  }

  const PartitionedHyperGraph& hypergraph_of_edge(const HyperedgeID e) const {
    int node = common::get_numa_node_of_edge(e);
    ASSERT(node < static_cast<int>(_hypergraphs.size()));
    return _hypergraphs[node];
  }

  PartitionedHyperGraph& hypergraph_of_edge(const HyperedgeID e) {
    return const_cast<PartitionedHyperGraph&>(static_cast<const NumaPartitionedHypergraph&>(*this).hypergraph_of_edge(e));
  }

  // ! Number of blocks
  PartitionID _k;
  // ! Hypergraph object around this partitioned hypergraph is wrapped
  Hypergraph* _hg;

  // ! Weight and size information for all blocks.
  parallel::scalable_vector<PartInfo> _part_info;
  // ! Partitioned NUMA Hypergraphs
  parallel::scalable_vector<PartitionedHyperGraph> _hypergraphs;
};

} // namespace ds
} // namespace mt_kahypar
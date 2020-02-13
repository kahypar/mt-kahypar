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

template <typename Hypergraph = Mandatory,
          bool track_border_vertices = true>
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

  class VertexPartInfo {
    public:
      explicit VertexPartInfo() :
        part_id(kInvalidPartition),
        num_incident_cut_hes(0UL) { }

    parallel::IntegralAtomicWrapper<PartitionID> part_id;
    parallel::IntegralAtomicWrapper<HypernodeID> num_incident_cut_hes;
  };

  using PinCountAtomic = parallel::IntegralAtomicWrapper<HypernodeID>;

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = false;
  static constexpr bool is_partitioned = true;

  explicit PartitionedHypergraph() :
    _k(0),
    _node(0),
    _hg(nullptr),
    _is_init_num_cut_hyperedges(false),
    _part_info(),
    _vertex_part_info(),
    _pins_in_part(),
    _connectivity_sets(0, 0) { }

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
    _k(k),
    _node(hypergraph.numaNode()),
    _hg(&hypergraph),
    _is_init_num_cut_hyperedges(false),
    _part_info(k),
    _vertex_part_info(hypergraph.initialNumNodes()),
    _pins_in_part(hypergraph.initialNumEdges() * k, PinCountAtomic(0)),
    _connectivity_sets(hypergraph.initialNumEdges(), k) { }

  explicit PartitionedHypergraph(const PartitionID k,
                                 const TaskGroupID,
                                 Hypergraph& hypergraph) :
    _k(k),
    _node(hypergraph.numaNode()),
    _hg(&hypergraph),
    _is_init_num_cut_hyperedges(false),
    _part_info(k),
    _vertex_part_info(),
    _pins_in_part(),
    _connectivity_sets(0, 0) {
    tbb::parallel_invoke([&] {
      _vertex_part_info.resize(hypergraph.initialNumNodes());
    }, [&] {
      _pins_in_part.assign(hypergraph.initialNumEdges() * k, PinCountAtomic(0));
    }, [&] {
      _connectivity_sets = ConnectivitySets(hypergraph.initialNumEdges(), k);
    });
  }

  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other) :
    _k(other._k),
    _node(other._node),
    _hg(other._hg),
    _is_init_num_cut_hyperedges(other._is_init_num_cut_hyperedges),
    _part_info(std::move(other._part_info)),
    _vertex_part_info(std::move(other._vertex_part_info)),
    _pins_in_part(std::move(other._pins_in_part)),
    _connectivity_sets(std::move(other._connectivity_sets)) { }

  PartitionedHypergraph & operator= (PartitionedHypergraph&& other) {
    _k = other._k;
    _node = other._node;
    _hg = other._hg;
    _is_init_num_cut_hyperedges = other._is_init_num_cut_hyperedges;
    _part_info = std::move(other._part_info);
    _vertex_part_info = std::move(other._vertex_part_info);
    _pins_in_part = std::move(other._pins_in_part);
    _connectivity_sets = std::move(other._connectivity_sets);
    return *this;
  }

  // ####################### General Hypergraph Stats #######################

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
    return _connectivity_sets.connectivitySet(local_id);
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
  bool setNodePart(const HypernodeID u, PartitionID id) {
    ASSERT(id != kInvalidPartition && id < _k);

    // Sets the node part of vertex u to id. The operation succeeds
    // if CAS operation on part_id of vertex u succeeds
    PartitionID invalid_id = kInvalidPartition;
    if ( vertexPartInfo(u).part_id.compare_and_exchange_strong(invalid_id, id) ) {

      // Update block weights
      ++_part_info[id].size;
      _part_info[id].weight += nodeWeight(u);

      // Update Pin Count Part of all incident edges
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        incrementPinCountInPart(he, id);
      }

      return true;
    } else {
      return false;
    }
  }

  // ! Sets the block id of an unassigned vertex u.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u,
                   PartitionID id,
                   parallel::scalable_vector<PartitionedHypergraph>& hypergraphs) {
    ASSERT(id != kInvalidPartition && id < _k);

    // Sets the node part of vertex u to id. The operation succeeds
    // if CAS operation on part_id of vertex u succeeds
    PartitionID invalid_id = kInvalidPartition;
    if ( vertexPartInfo(u).part_id.compare_and_exchange_strong(invalid_id, id) ) {

      // Update Pin Count Part of all incident edges
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        common::hypergraph_of_edge(he, hypergraphs).incrementPinCountInPart(he, id);
      }

      return true;
    } else {
      return false;
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  TRUE_SPECIALIZATION(track_border_vertices, bool)
  changeNodePart(const HypernodeID u,
                  PartitionID from,
                  PartitionID to,
                  const DeltaFunction& delta_func = NOOP_FUNC) {
    ASSERT(_hg->nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(from != kInvalidPartition && from < _k);
    ASSERT(to != kInvalidPartition && to < _k);
    ASSERT(from != to);
    ASSERT(_is_init_num_cut_hyperedges);

    if ( vertexPartInfo(u).part_id.compare_and_exchange_strong(from, to) ) {

      // Update block weights
      --_part_info[from].size;
      _part_info[from].weight -= nodeWeight(u);
      ++_part_info[to].size;
      _part_info[to].weight += nodeWeight(u);

      // Update Pin Count Part of all incident edges
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        HypernodeID pin_count_in_from_part_after = kInvalidHypernode;
        HypernodeID pin_count_in_to_part_after = kInvalidHypernode;
        updatePinCountOfHyperedge(he, *this, from, to, pin_count_in_from_part_after,
          pin_count_in_to_part_after, delta_func);

        bool no_pins_left_in_source_part = pin_count_in_from_part_after == 0;
        bool only_one_pin_in_to_part = pin_count_in_to_part_after == 1;

        const HypernodeID edge_size = edgeSize(he);
        if (no_pins_left_in_source_part && !only_one_pin_in_to_part &&
            pin_count_in_to_part_after == edge_size) {
          // In that case, hyperedge he becomes an internal hyperedge
          for (const HypernodeID& pin : pins(he)) {
            decrementIncidentNumCutHyperedges(pin);
          }
        } else if (!no_pins_left_in_source_part && only_one_pin_in_to_part &&
                  pin_count_in_from_part_after == edge_size - 1) {
          // In that case, hyperedge he becomes an cut hyperede
          for (const HypernodeID& pin : pins(he)) {
            incrementIncidentNumCutHyperedges(pin);
          }
        }
      }

      return true;
    } else {
      return false;
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  TRUE_SPECIALIZATION(track_border_vertices, bool)
  changeNodePart(const HypernodeID u,
                  PartitionID from,
                  PartitionID to,
                  parallel::scalable_vector<PartitionedHypergraph>& hypergraphs,
                  const DeltaFunction& delta_func = NOOP_FUNC) {
    ASSERT(_hg->nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(from != kInvalidPartition && from < _k);
    ASSERT(to != kInvalidPartition && to < _k);
    ASSERT(from != to);
    ASSERT(_is_init_num_cut_hyperedges);

    if ( vertexPartInfo(u).part_id.compare_and_exchange_strong(from, to) ) {

      // Update Pin Count Part of all incident edges
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        HypernodeID pin_count_in_from_part_after = kInvalidHypernode;
        HypernodeID pin_count_in_to_part_after = kInvalidHypernode;
        PartitionedHypergraph& hypergraph_of_he = common::hypergraph_of_edge(he, hypergraphs);
        updatePinCountOfHyperedge(he, hypergraph_of_he, from, to,
          pin_count_in_from_part_after, pin_count_in_to_part_after, delta_func);

        bool no_pins_left_in_source_part = pin_count_in_from_part_after == 0;
        bool only_one_pin_in_to_part = pin_count_in_to_part_after == 1;

        HypernodeID edge_size = hypergraph_of_he.edgeSize(he);
        if (no_pins_left_in_source_part && !only_one_pin_in_to_part &&
            pin_count_in_to_part_after == edge_size) {
          // In that case, hyperedge he becomes an internal hyperedge
          for (const HypernodeID& pin : hypergraph_of_he.pins(he)) {
            common::hypergraph_of_vertex(pin, hypergraphs).decrementIncidentNumCutHyperedges(pin);
          }
        } else if (!no_pins_left_in_source_part && only_one_pin_in_to_part &&
                  pin_count_in_from_part_after == edge_size - 1) {
          // In that case, hyperedge he becomes an cut hyperede
          for (const HypernodeID& pin : hypergraph_of_he.pins(he)) {
            common::hypergraph_of_vertex(pin, hypergraphs).incrementIncidentNumCutHyperedges(pin);
          }
        }
      }

      return true;
    } else {
      return false;
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  FALSE_SPECIALIZATION(track_border_vertices, bool)
  changeNodePart(const HypernodeID u,
                 PartitionID from,
                 PartitionID to,
                 const DeltaFunction& delta_func = NOOP_FUNC) {
    ASSERT(_hg->nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(from != kInvalidPartition && from < _k);
    ASSERT(to != kInvalidPartition && to < _k);
    ASSERT(from != to);
    ASSERT(_is_init_num_cut_hyperedges);

    if ( vertexPartInfo(u).part_id.compare_and_exchange_strong(from, to) ) {

      // Update block weights
      --_part_info[from].size;
      _part_info[from].weight -= nodeWeight(u);
      ++_part_info[to].size;
      _part_info[to].weight += nodeWeight(u);

      // Update Pin Count Part of all incident edges
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        updatePinCountOfHyperedge(he, *this, from, to, delta_func);
      }

      return true;
    } else {
      return false;
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  FALSE_SPECIALIZATION(track_border_vertices, bool)
  changeNodePart(const HypernodeID u,
                 PartitionID from,
                 PartitionID to,
                 parallel::scalable_vector<PartitionedHypergraph>& hypergraphs,
                 const DeltaFunction& delta_func = NOOP_FUNC) {
    ASSERT(_hg->nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(from != kInvalidPartition && from < _k);
    ASSERT(to != kInvalidPartition && to < _k);
    ASSERT(from != to);
    ASSERT(_is_init_num_cut_hyperedges);

    if ( vertexPartInfo(u).part_id.compare_and_exchange_strong(from, to) ) {

      // Update Pin Count Part of all incident edges
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        updatePinCountOfHyperedge(he, common::hypergraph_of_edge(he, hypergraphs),
          from, to, delta_func);
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
    return vertexPartInfo(u).part_id;
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  TRUE_SPECIALIZATION(track_border_vertices, bool) isBorderNode(const HypernodeID u) const {
    ASSERT(_is_init_num_cut_hyperedges);
    return vertexPartInfo(u).num_incident_cut_hes > 0;
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  FALSE_SPECIALIZATION(track_border_vertices, bool) isBorderNode(const HypernodeID u) const {
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( connectivity(he) > 1 ) {
        return true;
      }
    }
    return false;
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  FALSE_SPECIALIZATION(track_border_vertices, bool) isBorderNode(
    const HypernodeID u,
    const parallel::scalable_vector<PartitionedHypergraph>& hypergraphs) const {
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( common::hypergraph_of_edge(he, hypergraphs).connectivity(he) > 1 ) {
        return true;
      }
    }
    return false;
  }

  // ! Number of incident cut hyperedges of vertex u
  TRUE_SPECIALIZATION(track_border_vertices, HyperedgeID) numIncidentCutHyperedges(const HypernodeID u) const {
    return vertexPartInfo(u).num_incident_cut_hes;
  }

  // ! Number of incident cut hyperedges of vertex u
  FALSE_SPECIALIZATION(track_border_vertices, HyperedgeID) numIncidentCutHyperedges(const HypernodeID u) const {
    HyperedgeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      num_incident_cut_hyperedges += ( connectivity(he) > 1 );
    }
    return num_incident_cut_hyperedges;
  }

  // ! Number of incident cut hyperedges of vertex u
  FALSE_SPECIALIZATION(track_border_vertices, HyperedgeID) numIncidentCutHyperedges(
    const HypernodeID u,
    const parallel::scalable_vector<PartitionedHypergraph>& hypergraphs) const {
    HyperedgeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      num_incident_cut_hyperedges += ( common::hypergraph_of_edge(he, hypergraphs).connectivity(he) > 1 );
    }
    return num_incident_cut_hyperedges;
  }

  // ! Initializes the number of cut hyperedges for each vertex
  // ! NOTE, this function have to be called after initial partitioning
  // ! and before local search.
  TRUE_SPECIALIZATION(track_border_vertices, void) initializeNumCutHyperedges(
    const TaskGroupID task_group_id,
    const parallel::scalable_vector<PartitionedHypergraph>& hypergraphs = {}) {
    ASSERT(!_is_init_num_cut_hyperedges);
    _hg->doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
      ASSERT(partID(hn) != kInvalidPartition);
      for ( const HyperedgeID& he : incidentEdges(hn)) {
        const PartitionID he_connectivity = hypergraphs.empty() ? connectivity(he) :
          common::hypergraph_of_edge(he, hypergraphs).connectivity(he);
        if ( he_connectivity > 1 ) {
          incrementIncidentNumCutHyperedges(hn);
        }
      }
    });
    _is_init_num_cut_hyperedges = true;
  }

  // ! Initializes the number of cut hyperedges for each vertex
  // ! NOTE, this function have to be called after initial partitioning
  // ! and before local search.
  TRUE_SPECIALIZATION(track_border_vertices, void) initializeNumCutHyperedges(
    const parallel::scalable_vector<PartitionedHypergraph>& hypergraphs = {}) {
    ASSERT(!_is_init_num_cut_hyperedges);
    for ( const HypernodeID& hn : nodes() ) {
      ASSERT(partID(hn) != kInvalidPartition);
      for ( const HyperedgeID& he : incidentEdges(hn)) {
        const PartitionID he_connectivity = hypergraphs.empty() ? connectivity(he) :
          common::hypergraph_of_edge(he, hypergraphs).connectivity(he);
        if ( he_connectivity > 1 ) {
          incrementIncidentNumCutHyperedges(hn);
        }
      }
    }
    _is_init_num_cut_hyperedges = true;
  }

  // ! NOOP, in case border vertices should be not explicitly tracked
  FALSE_SPECIALIZATION(track_border_vertices, void) initializeNumCutHyperedges(
    const TaskGroupID,
    const parallel::scalable_vector<PartitionedHypergraph>& hypergraphs = {}) {
    unused(hypergraphs);
    _is_init_num_cut_hyperedges = true;
  }

  // ! NOOP, in case border vertices should be not explicitly tracked
  FALSE_SPECIALIZATION(track_border_vertices, void) initializeNumCutHyperedges(
    const parallel::scalable_vector<PartitionedHypergraph>& hypergraphs = {}) {
    unused(hypergraphs);
    _is_init_num_cut_hyperedges = true;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    return _connectivity_sets.connectivity(local_id);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID id) const {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(id != kInvalidPartition && id < _k);
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    return _pins_in_part[local_id * _k + id];;
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
    // Reset partition ids
    for ( VertexPartInfo& vertex_part_info : _vertex_part_info ) {
      vertex_part_info.part_id = kInvalidPartition;
      vertex_part_info.num_incident_cut_hes = 0;
    }

    // Reset pin count in part and connectivity set
    for ( const HyperedgeID& he : edges() ) {
      const HyperedgeID local_id = common::get_local_position_of_edge(he);
      for ( const PartitionID& block : connectivitySet(he) ) {
        _pins_in_part[static_cast<size_t>(local_id) * _k + block] = 0;
      }
      _connectivity_sets.clear(local_id);
    }

    // Reset block weights and sizes
    for ( PartitionID block = 0; block < _k; ++block ) {
      _part_info[block].weight = 0;
      _part_info[block].size = 0;
    }

    _is_init_num_cut_hyperedges = false;
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);

    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    _connectivity_sets.memoryConsumption(connectivity_set_node);

    parent->addChild("Part Info", sizeof(PartInfo) * _k);
    parent->addChild("Vertex Part Info", sizeof(VertexPartInfo) * _hg->initialNumNodes());
    parent->addChild("Pin Count In Part", sizeof(PinCountAtomic) * _k * _hg->initialNumEdges());
  }

 private:

  // ! Accessor for partition information of a vertex
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const VertexPartInfo& vertexPartInfo(const HypernodeID u) const {
    ASSERT(_hg->nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    HypernodeID local_id = common::get_local_position_of_vertex(u);
    ASSERT(common::get_numa_node_of_vertex(u) == _node,
           "Hypernode" << u << "is not part of numa node" << _node);
    ASSERT(local_id < _hg->initialNumNodes(), "Hypernode" << u << "does not exist");
    return _vertex_part_info[local_id];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE VertexPartInfo& vertexPartInfo(const HypernodeID u) {
    return const_cast<VertexPartInfo&>(static_cast<const PartitionedHypergraph&>(*this).vertexPartInfo(u));
  }

  // ####################### Partition Information #######################

  // ! Updates pin count in part if border vertices should be tracked.
  // ! The update process of the border vertices rely that
  // ! pin_count_in_from_part_after and pin_count_in_to_part_after are not reflecting
  // ! some intermediate state of the pin counts when several vertices move in parallel.
  // ! Therefore, we verify that the pin counts before the move sum up to the total edge size
  // ! and then verify that from part and to part are exactly decremented and incremented by one.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePinCountOfHyperedge(const HyperedgeID& he,
                                                                 PartitionedHypergraph& hypergraph_of_he,
                                                                 const PartitionID from,
                                                                 const PartitionID to,
                                                                 HypernodeID& pin_count_in_from_part_after,
                                                                 HypernodeID& pin_count_in_to_part_after,
                                                                 const DeltaFunction& delta_func) {
    HypernodeID pin_count_in_from_part_before = kInvalidHypernode;
    HypernodeID pin_count_in_to_part_before = kInvalidHypernode;
    pin_count_in_from_part_after = kInvalidHypernode;
    pin_count_in_to_part_after = kInvalidHypernode;
    HypernodeID edge_size = hypergraph_of_he.edgeSize(he);
    // In order to safely update the number of incident cut hyperedges and to compute
    // the delta of the move we need a stable snapshot of the pin count in from and to
    // part before and after the move. If we not do so, it can happen that due to concurrent
    // updates the pin count represents some intermediate state and the conditions
    // below are not triggered which leaves the data structure in an inconsistent
    // state. However, this should happen very rarely.
    while ( pin_count_in_from_part_before - 1 != pin_count_in_from_part_after ||
            pin_count_in_to_part_before + 1 != pin_count_in_to_part_after  ) {
      getPinCountInFromAndToPart(he, hypergraph_of_he, edge_size, from, to,
        pin_count_in_from_part_before, pin_count_in_to_part_before);
      pin_count_in_from_part_after =
        hypergraph_of_he.decrementPinCountInPart(he, from);
      pin_count_in_to_part_after =
        hypergraph_of_he.incrementPinCountInPart(he, to);

      if ( pin_count_in_from_part_before - 1 != pin_count_in_from_part_after ||
            pin_count_in_to_part_before + 1 != pin_count_in_to_part_after  ) {
        hypergraph_of_he.incrementPinCountInPart(he, from);
        hypergraph_of_he.decrementPinCountInPart(he, to);
        pin_count_in_from_part_after = kInvalidHypernode;
        pin_count_in_to_part_after = kInvalidHypernode;
      }
    }

    delta_func(he, hypergraph_of_he.edgeWeight(he), edge_size,
      pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  // ! Updates pin count in part if no border vertices should be tracked.
  // ! Note, in case delta_func is km1 or cut delta it can produce wrong deltas, because
  // ! pin_count_in_from_part_after and pin_count_in_to_part_after can reflect some
  // ! intermediate state of the pin counts when several vertices move in parallel.
  // ! To rely on the delta one should use the first update function.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePinCountOfHyperedge(const HyperedgeID& he,
                                                                 PartitionedHypergraph& hypergraph_of_he,
                                                                 const PartitionID from,
                                                                 const PartitionID to,
                                                                 const DeltaFunction& delta_func) {
    // If updated concurrently the pin counts in from and to part can represent some intermediate
    // state and it can happen that the below conditions are not triggered.
    HypernodeID pin_count_in_from_part_after = hypergraph_of_he.decrementPinCountInPart(he, from);
    HypernodeID pin_count_in_to_part_after = hypergraph_of_he.incrementPinCountInPart(he, to);

    // Note, this function can compute wrong results, see comment in
    // updatePinCountOfHyperedge(...)
    delta_func(he, hypergraph_of_he.edgeWeight(he), hypergraph_of_he.edgeSize(he),
      pin_count_in_from_part_after, pin_count_in_to_part_after);
  }

  // ! Decrements the number of incident cut hyperedges of hypernode u
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void decrementIncidentNumCutHyperedges(const HypernodeID u) {
    --vertexPartInfo(u).num_incident_cut_hes;
  }

  // ! Increments the number of incident cut hyperedges of hypernode u
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void incrementIncidentNumCutHyperedges(const HypernodeID u) {
    ++vertexPartInfo(u).num_incident_cut_hes;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID decrementPinCountInPart(const HyperedgeID e,
                                                                      const PartitionID id) {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(id != kInvalidPartition && id < _k);
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    const HypernodeID pin_count_after = --_pins_in_part[local_id * _k + id];
    if ( pin_count_after == 0 ) {
      // Connectivity of hyperedge decreased
      _connectivity_sets.remove(local_id, id);
    }
    return pin_count_after;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID incrementPinCountInPart(const HyperedgeID e,
                                                                      const PartitionID id) {
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(id != kInvalidPartition && id < _k);
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e),
           "Hyperedge" << e << "is not part of numa node" << _node);
    const HypernodeID pin_count_after = ++_pins_in_part[local_id * _k + id];
    if ( pin_count_after == 1 ) {
      // Connectivity of hyperedge increased
      _connectivity_sets.add(local_id, id);
    }
    return pin_count_after;
  }

  // ! Returns pin count in from and to part of hyperedge he such that
  // ! the pin counts represents a stable snapshot of the hyperedge
  // ! (required for safe change node part)
  void getPinCountInFromAndToPart(
    const HyperedgeID he,
    const PartitionedHypergraph& hypergraph_of_he,
    const HypernodeID edge_size,
    const PartitionID from,
    const PartitionID to,
    HypernodeID& pin_count_in_from_part_before,
    HypernodeID& pin_count_in_to_part_before) {
    HypernodeID sum_edge_size = 0;
    while ( sum_edge_size != edge_size ) {
      sum_edge_size = 0;
      pin_count_in_from_part_before = 0;
      pin_count_in_to_part_before = 0;
      for ( const PartitionID& block : hypergraph_of_he.connectivitySet(he) ) {
        const HypernodeID pin_count_in_part = hypergraph_of_he.pinCountInPart(he, block);
        sum_edge_size += pin_count_in_part;
        if ( block == from ) {
          pin_count_in_from_part_before = pin_count_in_part;
        } else if ( block == to ) {
          pin_count_in_to_part_before = pin_count_in_part;
        }
      }
    }
  }

  // ! Number of blocks
  PartitionID _k;
  // ! NUMA node of this partitioned hypergraph
  int _node;
  // ! Hypergraph object around this partitioned hypergraph is wrapped
  Hypergraph* _hg;

  // ! Indicates, if cut hyperedges are initialized
  bool _is_init_num_cut_hyperedges;
  // ! Weight and size information for all blocks.
  parallel::scalable_vector<PartInfo> _part_info;
  // ! Contains for each vertex the block and the number of
  // ! incident cut hyperedges
  parallel::scalable_vector<VertexPartInfo> _vertex_part_info;
  // ! For each hyperedge and each block, _pins_in_part stores the
  // ! number of pins in that block
  parallel::scalable_vector<PinCountAtomic> _pins_in_part;
  // ! For each hyperedge, _connectivity_sets stores the connectivity and the set of block ids
  // ! which the pins of the hyperedge belongs to
  ConnectivitySets _connectivity_sets;
};

} // namespace ds
} // namespace mt_kahypar
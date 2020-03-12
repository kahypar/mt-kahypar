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
#include "mt-kahypar/datastructures/partition_info.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory>
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
  using UnderlyingFactory = typename HypergraphFactory::UnderlyingFactory;
  using HardwareTopology = typename Hypergraph::HardwareTopology;
  using TBBNumaArena = typename Hypergraph::TBBNumaArena;
  using PartitionedHyperGraph = PartitionedHypergraph<UnderlyingHypergraph, UnderlyingFactory>;


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
  using BlockInfo = typename PartitionInfo::BlockInfo;

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = true;
  static constexpr bool is_partitioned = true;

  explicit NumaPartitionedHypergraph() { }

  explicit NumaPartitionedHypergraph(const PartitionID k,
                                     Hypergraph& hypergraph,
                                     vec<HypernodeWeight>& max_part_weight) :
    _k(k),
    _hg(&hypergraph),
    part_weight(k, 0),
    _hypergraphs() {
    for ( size_t node = 0; node < hypergraph.numNumaHypergraphs(); ++node) {
      _hypergraphs.emplace_back(k, hypergraph.numaHypergraph(node));
    }
  }

  explicit NumaPartitionedHypergraph(const PartitionID k,
                                     const TaskGroupID task_group_id,
                                     Hypergraph& hypergraph,
                                     vec<HypernodeWeight>& max_part_weight) :
    _k(k),
    _hg(&hypergraph),
    part_weight(k, 0),
    _hypergraphs()
  {
    TBBNumaArena::instance().execute_sequential_on_all_numa_nodes(task_group_id, [&](const int) {
      _hypergraphs.emplace_back();
    });

    // Construct NUMA partitioned hypergraphs in parallel
    TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
      _hypergraphs[node] = PartitionedHyperGraph(k, task_group_id, hypergraph.numaHypergraph(node), max_part_weight);
    });
  }

  NumaPartitionedHypergraph(const NumaPartitionedHypergraph&) = delete;
  NumaPartitionedHypergraph & operator= (const NumaPartitionedHypergraph &) = delete;

  NumaPartitionedHypergraph(NumaPartitionedHypergraph&& other) = default;

  NumaPartitionedHypergraph & operator= (NumaPartitionedHypergraph&& other) = default;

  ~NumaPartitionedHypergraph() {
    freeInternalData();
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
    PartitionedHyperGraph& hypergraph_of_he = hypergraph_of_edge(he);
    for ( const HypernodeID& pin : pins(he) ) {
      hypergraph_of_he.incrementPinCountInPartWithoutGainUpdate(he, partID(pin));
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

  // ####################### Partition Information #######################

  /**
   *
   * No two threads are allowed to move the same vertex concurrently.
   *
   * The functions setNodePart(...) and changeNodePart(...) are implemented thread-safe
   * and lock-free. We update the following information when moving a node:
   *  1.) Block id of a vertex: part[ .. ] / partID(.. )
   *  2.) Pin count in a block of a hyperedge: pins_in_part[he * k + block] / pinCountInPart(he, block)
   *  3.) Connectivity set of a hyperedge (spanned blocks): connectivitySet(e, block), connectivity(he)
   *  4.) Block weights. part_weight[ .. ] / partWeight( .. )
   *  5.) The sum of incident edge weights of a node with zero/one pin in a block: moveFromBenefit(u, block), moveToPenalty(u, block)
   *
   * 2.) - 5.) are based on atomics.
   * In a parallel setting, a decently accurate km1 gain can be computed as moveFromBenefit(u, partID(u)) - moveToPenalty(u, target_block)
   * 2.) is used to maintain 5.)
   * 3.) is a commodity used for initial partitioning
   *
   * There is a version of setNodePart and changeNodePart that checks whether the move is feasible with
   * respect to the balance constraint, by trying to subtract the node weight from a user-supplied budget
   * for the destination block. If that budget falls below zero, the move fails and the budget update is reverted.
   * Otherwise the move succeeds and we update the remaining information.
   * The idea behind maintaining a budget is to alleviate contention on part weights.
   * Instead, each thread can maintain its own budget and 'steal' some from other threads if it needs to.
   * Therefore, changeNodePartWithBalanceCheckAndGainUpdates does not maintain part weights.
   * They can easily be recovered by
   *    a) thread-local part weight deltas and calling applyPartWeightUpdates(..)
   *    or
   *    b) the overall budget for a block across all threads is the max part weight minus the current part weight

   * It is guaranteed, that if all concurrent vertex move operations are finished, the partition information
   * is in a correct state. The only exception is the above mentioned postprocessing necessary for changeNodePartWithBalanceCheck.
   *
   */

  // ! Sets the block id of an unassigned vertex u.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  bool setNodePart(const HypernodeID u, const PartitionID id) {
    if ( hypergraph_of_vertex(u).setNodePart(u, id, _hypergraphs) ) {
      // Update block weights
      part_weight[id] += nodeWeight(u);
      return true;
    } else {
      return false;
    }
  }

  // ! Sets the block id of an unassigned vertex u.
  // ! Returns true, if assignment of unassigned vertex u to block id succeeds.
  // ! Note, that in contrast to setNodePart the block weights and sizes and also
  // ! the pin count in part of all incident hyperedges is not updated. In order to
  // ! update those stats, one has to call initializePartition(...) after all
  // ! block ids are assigned.
  bool setOnlyNodePart(const HypernodeID u, PartitionID id) {
    return hypergraph_of_vertex(u).setOnlyNodePart(u, id);
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    if ( hypergraph_of_vertex(u).changeNodePart(u, from, to, _hypergraphs, delta_func) ) {
      // Update block weights
      part_weight[from] -= nodeWeight(u);
      part_weight[to] += nodeWeight(u);
      return true;
    } else {
      return false;
    }
  }

  // TODO at the moment there is no implementation of changeNodePartWithBalanceCheckAndGainUpdate( .. )
  // and the corresponding constant time gain queries
  // because the FM implementation is not yet ready for the numa-aware mode
  // this will be added once FM works properly in non-numa-aware mode

  // ! Block which vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    return hypergraph_of_vertex(u).partID(u);
  }


  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, part info, pin counts in part and border
  // ! vertices have to be computed in a postprocessing step.
  void initializePartition(const TaskGroupID task_group_id) {
    tbb::parallel_invoke([&] {
      // Compute Part Info
      TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int socket) {
        _hypergraphs[socket].initializeBlockWeights();
      });
      for (const auto& socket_phg : _hypergraphs) {
        for (PartitionID p = 0; p < _k; ++p) {
          part_weight[p] += socket_phg.partWeight(p);
        }
      }
    }, [&] {
      // Compute Pin Counts
      TBBNumaArena::instance().execute_parallel_on_all_numa_nodes(task_group_id, [&](const int node) {
        _hypergraphs[node].initializePartition(task_group_id, _hypergraphs);
      });
    });
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    return hypergraph_of_vertex(u).isBorderNode(u);
  }

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    return hypergraph_of_vertex(u).numIncidentCutHyperedges(u);
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
    return part_weight[id].load(std::memory_order_relaxed);
  }


  // ! Reset partition (not thread-safe)
  void resetPartition() {
    for ( PartitionedHyperGraph& partitioned_hypergraph : _hypergraphs  ) {
      partitioned_hypergraph.resetPartition();
    }

    // Reset block weights and sizes
    for ( PartitionID block = 0; block < _k; ++block ) {
      part_weight[block].store(0);
    }
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      tbb::parallel_for(0UL, _hypergraphs.size(), [&](const size_t i) {
        _hypergraphs[i].freeInternalData();
      }, tbb::static_partitioner());
    }
    _k = 0;
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Part Info", sizeof(BlockInfo) * _k);
    for ( size_t node = 0; node < _hypergraphs.size(); ++node ) {
      utils::MemoryTreeNode* numa_hypergraph_node = parent->addChild(
        "NUMA Partitioned Hypergraph " + std::to_string(node));
      _hypergraphs[node].memoryConsumption(numa_hypergraph_node);
    }
  }

  // ####################### Extract Block #######################

  // ! Extracts a block of a partition as separate hypergraph.
  // ! It also returns an mapping, that maps each vertex of the original
  // ! hypergraph to a vertex in the extracted hypergraph.
  // ! Furthermore, if cut_net_splitting is activated, hyperedges that
  // ! contains more than one block are splitted and otherwise only
  // ! only nets that are internal in 'block' are in the extracted hypergraph
  std::pair<Hypergraph, parallel::scalable_vector<HypernodeID> > extract(
    const TaskGroupID& task_group_id,
    const PartitionID block,
    const bool cut_net_splitting) {
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
    parallel::scalable_vector<int> vertices_to_numa_node;
    parallel::scalable_vector<HyperedgeWeight> hyperedge_weight;
    parallel::scalable_vector<HypernodeWeight> hypernode_weight;
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weight.resize(num_hyperedges);
      doParallelForAllEdges(task_group_id, [&](const HyperedgeID he) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          const HyperedgeID extracted_edge_id = he_mapping[originalEdgeID(he)];
          hyperedge_weight[extracted_edge_id] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == block ) {
              edge_vector[extracted_edge_id].push_back(hn_mapping[originalNodeID(pin)]);
            }
          }
        }
      });
    }, [&] {
      hypernode_weight.resize(num_hypernodes);
      vertices_to_numa_node.resize(num_hypernodes);
      doParallelForAllNodes(task_group_id, [&](const HypernodeID hn) {
        if ( partID(hn) == block ) {
          const HypernodeID extracted_vertex_id = hn_mapping[originalNodeID(hn)];
          hypernode_weight[extracted_vertex_id] = nodeWeight(hn);
          vertices_to_numa_node[extracted_vertex_id] = common::get_numa_node_of_vertex(hn);
        }
      });
    });

    // Construct hypergraph
    Hypergraph extracted_hypergraph = HypergraphFactory::construct(
      task_group_id, num_hypernodes, num_hyperedges,
      edge_vector, std::move(vertices_to_numa_node),
      hyperedge_weight.data(), hypernode_weight.data());

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
  PartitionID _k = 0;
  // ! Hypergraph object around this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Weight and information for all blocks.
  vec< CAtomic<HypernodeWeight> > part_weight;

  // ! Partitioned NUMA Hypergraphs
  parallel::scalable_vector<PartitionedHyperGraph> _hypergraphs;
};

} // namespace ds
} // namespace mt_kahypar
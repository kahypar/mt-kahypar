/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

  // REVIEW NOTE do we still need these templates. we only have one hypergraph type now, right?

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory>
class PartitionedHypergraph {
private:
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
  using CommunityIterator = typename Hypergraph::CommunityIterator;

  using AtomicFlag = parallel::IntegralAtomicWrapper<bool>;
  template<typename T>
  using ThreadLocalVector = tbb::enumerable_thread_specific<parallel::scalable_vector<T>>;

  // ! Function that will be called for each incident hyperedge of a moved vertex with the following arguments
  // !  1) hyperedge ID, 2) weight, 3) size, 4) pin count in from-block after move, 5) pin count in to-block after move
  // ! Can be implemented to obtain correct km1 or cut improvements of the move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  // REVIEW NOTE: Can't we use a lambda in changeNodePart. And write a second function that calls the first with a lambda that does nothing.
  // Then we could guarantee inlining
  // This would also reduce the code/documentation copy-pasta for with or without gain updates

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_partitioned = true;

  PartitionedHypergraph() = default;

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
    _k(k),
    _hg(&hypergraph),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(
        "Refinement", "part_ids", hypergraph.initialNumNodes(), false, false),
    _pins_in_part(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize(), false),
    _connectivity_set(hypergraph.initialNumEdges(), k, false),
    _move_to_penalty(
        "Refinement", "move_to_penalty", hypergraph.initialNumNodes() * k, true, false),
    _move_from_benefit(
        "Refinement", "move_from_benefit", hypergraph.initialNumNodes(), true, false),
    _pin_count_update_ownership(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true, false) {
    _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition, false);
  }

  explicit PartitionedHypergraph(const PartitionID k,
                                 const TaskGroupID,
                                 Hypergraph& hypergraph) :
    _k(k),
    _hg(&hypergraph),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(),
    _pins_in_part(),
    _connectivity_set(0, 0),
    _move_to_penalty(),
    _move_from_benefit(),
    _pin_count_update_ownership() {
    tbb::parallel_invoke([&] {
      _part_ids.resize(
        "Refinement", "vertex_part_info", hypergraph.initialNumNodes());
      _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition);
    }, [&] {
      _pins_in_part.initialize(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize());
    }, [&] {
      _connectivity_set = ConnectivitySets(hypergraph.initialNumEdges(), k);
    }, [&] {
      _move_to_penalty.resize(
        "Refinement", "move_to_penalty", hypergraph.initialNumNodes() * k, true);
    }, [&] {
      _move_from_benefit.resize(
        "Refinement", "move_from_benefit", hypergraph.initialNumNodes(), true);
    }, [&] {
      _pin_count_update_ownership.resize(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true);
    });
  }

  // REVIEW NOTE why do we delete copy assignment/construction? wouldn't it be useful to make a copy, e.g. for initial partitioning
  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other) = default;
  PartitionedHypergraph & operator= (PartitionedHypergraph&& other) = default;

  ~PartitionedHypergraph() {
    freeInternalData();
  }

  // ####################### General Hypergraph Stats ######################

  Hypergraph& hypergraph() {
    ASSERT(_hg);
    return *_hg;
  }

  void setHypergraph(Hypergraph& hypergraph) {
    _hg = &hypergraph;
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _hg->initialNumNodes();
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _hg->numRemovedHypernodes();
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _hg->initialNumEdges();
  }

  HyperedgeID numGraphEdges() const {
    return _hg->numGraphEdges();
  }

  HyperedgeID numNonGraphEdges() const {
    return _hg->numNonGraphEdges();
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _hg->initialNumPins();
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _hg->initialTotalVertexDegree();
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
  void doParallelForAllNodes(const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    _hg->doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllEdges(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    _hg->doParallelForAllEdges(f);
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return _hg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
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
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    return _connectivity_set.connectivitySet(e);
  }

  // ####################### Hypernode Information #######################

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
      _part_weights[block] += delta;
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

  HyperedgeID graphEdgeID(const HyperedgeID e) const {
    return _hg->graphEdgeID(e);
  }

  HyperedgeID nonGraphEdgeID(const HyperedgeID e) const {
    return _hg->nonGraphEdgeID(e);
  }

  HypernodeID graphEdgeHead(const HyperedgeID e, const HypernodeID tail) const {
    return _hg->graphEdgeHead(e, tail);
  }

  // ####################### Partition Information #######################

  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(_part_ids[u] == kInvalidPartition);
    _part_ids[u] = p;
  }

  void setNodePart(const HypernodeID u, PartitionID p) {
    setOnlyNodePart(u, p);
    _part_weights[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
    for (HyperedgeID he : incidentEdges(u)) {
      incrementPinCountInPartWithoutGainUpdate(he, p);
    }
  }

  // This function type does not update part weights but instead updates a slack
  // If this is used, part weights have to be recomputed at the end of partition initialization
  // This function serves only as the counterpart to the more useful changeNodePartWithBalanceCheck function
  // The use of thread-local slack variables is supposed to alleviate contention on the part_weight vector
  bool setOnlyNodePartWithBalanceCheck(const HypernodeID u, PartitionID p, CAtomic<HypernodeWeight>& budget_p) {
    const HypernodeWeight wu = nodeWeight(u);
    if (budget_p.sub_fetch(wu, std::memory_order_relaxed) >= 0) {
      _part_ids[u] = p;
      return true;
    } else {
      budget_p.fetch_add(wu, std::memory_order_relaxed);
      return false;
    }
  }

  void changeOnlyNodePart(const HypernodeID u, PartitionID from, PartitionID to) {
    nodeGainAssertions(u, from);
    nodeGainAssertions(u, to);
    ASSERT(from != to);

    const HypernodeWeight wu = nodeWeight(u);
    _part_weights[to].fetch_add(wu, std::memory_order_relaxed);
    _part_weights[from].fetch_sub(wu, std::memory_order_relaxed);
    _part_ids[u] = to;
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u, PartitionID from, PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    changeOnlyNodePart(u, from,  to);
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      // REVIEW NOTE wouldn't it be more elegant to write this with a spinlock directly?
      while ( !updatePinCountOfHyperedgeWithoutGainUpdates(he, from, to, delta_func) );
    }
    return true;
  }

  template<typename F>
  bool changeNodePartFullUpdate(
          const HypernodeID u, PartitionID from, PartitionID to, HypernodeWeight max_weight_to, F&& report_success) {
    assert(partID(u) == from);
    assert(from != to);
    const HypernodeWeight wu = nodeWeight(u);
    const HypernodeWeight to_weight_after = _part_weights[to].add_fetch(wu, std::memory_order_relaxed);
    const HypernodeWeight from_weight_after = _part_weights[from].fetch_sub(wu, std::memory_order_relaxed);
    if (to_weight_after <= max_weight_to & from_weight_after > 0) {
      report_success();
      _part_ids[u] = to;
      for (HyperedgeID he: incidentEdges(u)) {
        while ( !updatePinCountOfHyperedgeWithGainUpdates(he, from, to) );
      }
      return true;
    } else {
      _part_weights[to].fetch_sub(wu, std::memory_order_relaxed);
      _part_weights[from].fetch_add(wu, std::memory_order_relaxed);
      return false;
    }
  }


  // Additionally rejects the requested move if it violates balance
  // must recompute part_weights once finished moving nodes
  bool changeNodePartWithBalanceCheckAndGainUpdatesWithoutPartWeightUpdates(const HypernodeID u,
                                                                            PartitionID from,
                                                                            CAtomic<HypernodeWeight>& budget_from,
                                                                            PartitionID to,
                                                                            CAtomic<HypernodeWeight>& budget_to) {
    const bool success = setOnlyNodePartWithBalanceCheck(u, to, budget_to);
    if (success) {
      budget_from.fetch_add(nodeWeight(u), std::memory_order_relaxed);
      for (HyperedgeID he : incidentEdges(u)) {
        while ( !updatePinCountOfHyperedgeWithGainUpdates(he, from, to) );
      }
    }
    return success;
  }

  void setPartWeight(PartitionID p, HypernodeWeight w) {
    _part_weights[p].store(w);
  }

  // ! Only for testing
  void recomputePartWeights() {
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].store(0);
    }

    for (HypernodeID u : nodes()) {
      _part_weights[ partID(u) ] += nodeWeight(u);
    }
  }

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    return _part_ids[u];
  }


  // requires pinCountInPart to be computed
  void initializeGainInformation() {
    // check whether part has been initialized
    ASSERT([&] {
      if (_part_ids.size() != initialNumNodes()) {
        return false;
      }
      for (HypernodeID u : nodes()) {
        if (partID(u) == kInvalidPartition || partID(u) > k()) {
          return false;
        }
      }
      return true;
    } ());

    // we can either assume that pinCountInPart has been initialized and use this information
    // or recompute the gain information from scratch, which would be independent but requires more memory volume
    // since the LP refiner currently does not use and update the gain information, we rely on pinCountInPart being up to date
    // as we have to call initializeGainInformation() before every FM call
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_mtp(_k, 0);

    auto accumulate_and_assign = [&](tbb::blocked_range<HypernodeID>& r) {

      vec<HyperedgeWeight>& l_move_to_penalty = ets_mtp.local();

      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( nodeIsEnabled(u) ) {
          const PartitionID from = partID(u);
          HyperedgeWeight incident_edges_weight = 0;
          HyperedgeWeight l_move_from_benefit = 0;
          for (HyperedgeID he : incidentEdges(u)) {

            HyperedgeWeight edge_weight = edgeWeight(he);
            if (pinCountInPart(he, from) == 1) {
              l_move_from_benefit += edge_weight;
            }

            for (const PartitionID block : connectivitySet(he)) {
              l_move_to_penalty[block] -= edge_weight;
            }
            incident_edges_weight += edge_weight;
          }

          _move_from_benefit[u].store(l_move_from_benefit, std::memory_order_relaxed);
          for (PartitionID p = 0; p < _k; ++p) {
            _move_to_penalty[u * _k + p].store(l_move_to_penalty[p] + incident_edges_weight, std::memory_order_relaxed);
            l_move_to_penalty[p] = 0;
          }
        }
      }

    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()), accumulate_and_assign);
  }

  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, part info and pin counts in part
  void initializePartition(const TaskGroupID ) {
    tbb::parallel_invoke(
            [&] { initializeBlockWeights(); },
            [&] { initializePinCountInPart(); }
    );
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

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    HypernodeID incident_cut_hes = 0;
    for (const HyperedgeID he : incidentEdges(u)) {
      incident_cut_hes += connectivity(he) > 1 ? 1 : 0;
    }
    return incident_cut_hes;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    return _connectivity_set.connectivity(e);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    return _pins_in_part.pinCountInPart(e, p);
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_weights[p].load(std::memory_order_relaxed);
  }

  HyperedgeWeight moveFromBenefit(const HypernodeID u) const {
    return _move_from_benefit[u].load(std::memory_order_relaxed);
  }

  HyperedgeWeight moveToPenalty(const HypernodeID u, PartitionID p) const {
    return _move_to_penalty[u * _k + p].load(std::memory_order_relaxed);
  }

  HyperedgeWeight km1Gain(const HypernodeID u, PartitionID from, PartitionID to) const {
    unused(from);
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveFromBenefit(u) - moveToPenalty(u, to);
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    _part_ids.assign(_part_ids.size(), kInvalidPartition, false);
    for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);

    // Reset pin count in part and connectivity set
    for ( const HyperedgeID& he : edges() ) {
      for ( const PartitionID& block : connectivitySet(he) ) {
        _pins_in_part.setPinCountInPart(he, block, 0);
      }
      _connectivity_set.clear(he);
    }
  }

  // ! Only for testing
  HyperedgeWeight moveFromBenefitRecomputed(const HypernodeID u) const {
    const PartitionID p = partID(u);
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (pinCountInPart(e, p) == 1) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  // ! Only for testing
  HyperedgeWeight moveToPenaltyRecomputed(const HypernodeID u, PartitionID p) const {
    HyperedgeWeight w = 0;
    for (HyperedgeID e : incidentEdges(u)) {
      if (pinCountInPart(e, p) == 0) {
        w += edgeWeight(e);
      }
    }
    return w;
  }

  // ! Only for testing
  void recomputeMoveFromBenefit(const HypernodeID u) {
    _move_from_benefit[u].store(moveFromBenefitRecomputed(u));
  }

  // ! Only for testing
  bool checkTrackedPartitionInformation() {
    for (HyperedgeID e : edges()) {
      for (PartitionID i = 0; i < k(); ++i) {
        assert(pinCountInPart(e, i) == pinCountInPartRecomputed(e, i));
      }
    }
    for (HypernodeID u : nodes()) {
      assert(moveFromBenefit(u) == moveFromBenefitRecomputed(u));
      for (PartitionID i = 0; i < k(); ++i) {
        if (partID(u) != i) {
          assert(moveToPenalty(u, i) == moveToPenaltyRecomputed(u, i));
        }
      }
    }
    return true;
  }

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
// TODO finish this function when everything else is done

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);
    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    _connectivity_set.memoryConsumption(connectivity_set_node);

    parent->addChild("Part Info", sizeof(CAtomic<HypernodeWeight>) * _k);
    parent->addChild("Vertex Part Info", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("Pin Count In Part", _pins_in_part.size_in_bytes());
    parent->addChild("HE Ownership", sizeof(AtomicFlag) * _hg->initialNumNodes());
  }


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
          hn_mapping[hn] = num_hypernodes++;
        }
      }
    }, [&] {
      for ( const HyperedgeID& he : edges() ) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          he_mapping[he] = num_hyperedges++;
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
      doParallelForAllEdges([&](const HyperedgeID he) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          hyperedge_weight[he_mapping[he]] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == block ) {
              edge_vector[he_mapping[he]].push_back(hn_mapping[pin]);
            }
          }
        }
      });
    }, [&] {
      hypernode_weight.resize(num_hypernodes);
      doParallelForAllNodes([&](const HypernodeID hn) {
        if ( partID(hn) == block ) {
          hypernode_weight[hn_mapping[hn]] = nodeWeight(hn);
        }
      });
    });

    // Construct hypergraph
    Hypergraph extracted_hypergraph = HypergraphFactory::construct(
            task_group_id, num_hypernodes, num_hyperedges,
            edge_vector, hyperedge_weight.data(), hypernode_weight.data());

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( partID(hn) == block ) {
        const HypernodeID extracted_hn = hn_mapping[hn];
        extracted_hypergraph.setCommunityID(extracted_hn, _hg->communityID(hn));
      }
    });
    extracted_hypergraph.initializeCommunities();

    return std::make_pair(std::move(extracted_hypergraph), std::move(hn_mapping));
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      tbb::parallel_invoke( [&] {
        parallel::parallel_free(_part_ids, _pin_count_update_ownership);
      }, [&] {
        parallel::free(_pins_in_part.data());
      }, [&] {
        _connectivity_set.freeInternalData();
      } );
    }
    _k = 0;
  }

 private:

  void applyPartWeightUpdates(vec<HypernodeWeight>& part_weight_deltas) {
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].fetch_add(part_weight_deltas[p], std::memory_order_relaxed);
    }
  }

  void initializeBlockWeights() {
    auto accumulate = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HypernodeWeight> pws(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( nodeIsEnabled(u) ) {
          const PartitionID pu = partID( u );
          const HypernodeWeight wu = nodeWeight( u );
          pws[pu] += wu;
        }
      }
      applyPartWeightUpdates(pws);
    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  void initializePinCountInPart() {
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);

    auto assign = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HypernodeID>& pin_counts = ets_pin_count_in_part.local();
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {

        for (HypernodeID pin : pins(he)) {
          ++pin_counts[partID(pin)];
        }

        for (PartitionID p = 0; p < _k; ++p) {
          assert(pinCountInPart(he, p) == 0);
          if (pin_counts[p] > 0) {
            _connectivity_set.add(he, p);
            _pins_in_part.setPinCountInPart(he, p, pin_counts[p]);
          }
          pin_counts[p] = 0;
        }

      }
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HyperedgeID(0), initialNumEdges()), assign);
  }

  HypernodeID pinCountInPartRecomputed(const HyperedgeID e, PartitionID p) const {
    HypernodeID pcip = 0;
    for (HypernodeID u : pins(e)) {
      if (partID(u) == p) {
        pcip++;
      }
    }
    return pcip;
  }

  void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
    unused(u);
    unused(p);
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(u * _k + p < _move_to_penalty.size());
    ASSERT(u < _move_from_benefit.size());
  }

  // ! Updates pin count in part if border vertices should be tracked.
  // ! The update process of the border vertices rely that
  // ! pin_count_in_from_part_after and pin_count_in_to_part_after are not reflecting
  // ! some intermediate state of the pin counts when several vertices move in parallel.
  // ! Therefore, the current thread, which tries to modify the pin counts of the hyperedge,
  // ! try to acquire the ownership of the hyperedge and on success, pin counts are updated.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool updatePinCountOfHyperedgeWithoutGainUpdates(const HyperedgeID& he,
                                                                                   const PartitionID from,
                                                                                   const PartitionID to,
                                                                                   const DeltaFunction& delta_func) {
    // In order to safely update the number of incident cut hyperedges and to compute
    // the delta of a move we need a stable snapshot of the pin count in from and to
    // part before and after the move. If we not do so, it can happen that due to concurrent
    // updates the pin count represents some intermediate state and the conditions
    // below are not triggered which leaves the data structure in an inconsistent
    // state. However, this should happen very rarely.
    bool expected = 0;
    bool desired = 1;
    ASSERT(he < _pin_count_update_ownership.size());
    if ( _pin_count_update_ownership[he].compare_exchange_strong(expected, desired, std::memory_order_acq_rel) ) {
      // In that case, the current thread acquires the ownership of the hyperedge and can
      // safely update the pin counts in from and to part.
      const HypernodeID pin_count_in_from_part_after = decrementPinCountInPartWithoutGainUpdate(he, from);
      const HypernodeID pin_count_in_to_part_after = incrementPinCountInPartWithoutGainUpdate(he, to);
      delta_func(he, edgeWeight(he), edgeSize(he),
        pin_count_in_from_part_after, pin_count_in_to_part_after);
      _pin_count_update_ownership[he].store(false, std::memory_order_acq_rel);
      return true;
    }

    return false;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _pins_in_part.decrementPinCountInPart(e, p);
    if ( pin_count_after == 0 ) {
      _connectivity_set.remove(e, p);
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _pins_in_part.incrementPinCountInPart(e, p);
    if ( pin_count_after == 1 ) {
      _connectivity_set.add(e, p);
    }
    return pin_count_after;
  }


  // REVIEW NOTE documentation duplicated. no more copy-pasta please
  // if you can guarantee that the function call for delta_func is inlined, it's even the same code

  // ! Updates pin count in part if border vertices should be tracked.
  // ! The update process of the border vertices rely that
  // ! pin_count_in_from_part_after and pin_count_in_to_part_after are not reflecting
  // ! some intermediate state of the pin counts when several vertices move in parallel.
  // ! Therefore, the current thread, which tries to modify the pin counts of the hyperedge,
  // ! try to acquire the ownership of the hyperedge and on success, pin counts are updated.
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool updatePinCountOfHyperedgeWithGainUpdates(const HyperedgeID& he,
                                                                                const PartitionID from,
                                                                                const PartitionID to) {
    // In order to safely update the number of incident cut hyperedges and to compute
    // the delta of a move we need a stable snapshot of the pin count in from and to
    // part before and after the move. If we not do so, it can happen that due to concurrent
    // updates the pin count represents some intermediate state and the conditions
    // below are not triggered which leaves the data structure in an inconsistent
    // state. However, this should happen very rarely.
    bool expected = 0;
    bool desired = 1;
    ASSERT(he < _pin_count_update_ownership.size());
    if ( _pin_count_update_ownership[he].compare_exchange_strong(expected, desired, std::memory_order_acq_rel) ) {
      // In that case, the current thread acquires the ownership of the hyperedge and can
      // safely update the pin counts in from and to part.
      decrementPinCountInPartWithGainUpdate(he, from);
      incrementPinCountInPartWithGainUpdate(he, to);
      _pin_count_update_ownership[he].store(false, std::memory_order_acq_rel);
      return true;
    }

    return false;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPartWithGainUpdate(const HyperedgeID e, const PartitionID p) {
    const HypernodeID pin_count_after = decrementPinCountInPartWithoutGainUpdate(e, p);
    if (pin_count_after == 1) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        if (partID(u) == p)
          _move_from_benefit[u].fetch_add(we, std::memory_order_relaxed);
      }
    } else if (pin_count_after == 0) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        _move_to_penalty[u *_k + p].fetch_add(we, std::memory_order_relaxed);
      }
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPartWithGainUpdate(const HyperedgeID e, const PartitionID p) {
    const HypernodeID pin_count_after = incrementPinCountInPartWithoutGainUpdate(e, p);
    if (pin_count_after == 1) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        _move_to_penalty[u * _k + p].fetch_sub(we, std::memory_order_relaxed);
      }
    } else if (pin_count_after == 2) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        if (partID(u) == p)
          _move_from_benefit[u].fetch_sub(we, std::memory_order_relaxed);
      }
    }
    return pin_count_after;
  }

  // ! Number of blocks
  PartitionID _k = 0;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Weight and information for all blocks.
  vec< CAtomic<HypernodeWeight> > _part_weights;

  // ! Current block IDs of the vertices
  Array< PartitionID > _part_ids;

  // ! For each hyperedge and each block, _pins_in_part stores the
  // ! number of pins in that block
  PinCountInPart _pins_in_part;

  // TODO we probably don't need connectivity sets any more, except in the IP hypergraphs which don't need parallelism support
  // ! For each hyperedge, _connectivity_set stores the set of blocks that the hyperedge spans
  ConnectivitySets _connectivity_set;

  // ! For each node and block, the sum of incident edge weights with zero pins in that part
  Array< CAtomic<HyperedgeWeight> > _move_to_penalty;

  // ! For each node and block, the sum of incident edge weights with exactly one pin in that part
  Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

  // ! In order to update the pin count of a hyperedge thread-safe, a thread must acquire
  // ! the ownership of a hyperedge via a CAS operation.
  Array<AtomicFlag> _pin_count_update_ownership;
};

} // namespace ds
} // namespace mt_kahypar
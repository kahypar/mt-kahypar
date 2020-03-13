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
#include <external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h>

#include "tbb/parallel_invoke.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {


/**
 * A note on ID mappings for nodes and hyperedges.
 * There are three types of IDs: 1) global IDs 2) local IDs 3) original IDs
 *
 * Original IDs are consecutive in [0, |V|) or [0, |E|)
 * Local IDs are consecutive for the stuff allocated on each NUMA node/socket. These should be used to access arrays
 * in this class.
 * Global IDs consist of the local IDs in the bottom 48 bits and the NUMA/socket ID in the top 16 bits.
 *
 * All functions are supposed to handle global IDs, which is why they must perform a call to
 * common::get_local_position_of_{vertex/edge} before accessing memory.
 * All functions on this hypergraph can also be called with local IDs as the global-to-local mapping only hides
 * the NUMA ID bits.
 * TODO: this is currently not true! the assertions assume that they're given global IDs
 *
 * The pin and incident hyperedge IDs stored in _hg are global IDs, so beware before using them to access arrays!
 *
 */

// Forward
template <typename Hypergraph,
          typename HypergraphFactory>
class NumaPartitionedHypergraph;

template <typename Hypergraph = Mandatory,
          typename HypergraphFactory = Mandatory>
class PartitionedHypergraph {
private:
  template <typename HyperGraph,
          typename HyperGraphFactory>
  friend class NumaPartitionedHypergraph;
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");
  static_assert(!Hypergraph::is_numa_aware,  "Only non-numa-aware hypergraphs are allowed");

  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
  using CommunityIterator = typename Hypergraph::CommunityIterator;

  // ! Function that will be called for each incident hyperedge of a moved vertex with the following arguments
  // !  1) hyperedge ID, 2) weight, 3) size, 4) pin count in from-block after move, 5) pin count in to-block after move
  // ! Can be implemented to obtain correct km1 or cut improvements of the move
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }


 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_numa_aware = false;
  static constexpr bool is_partitioned = true;

  explicit PartitionedHypergraph() : connectivity_sets(0, 0) { }

  explicit PartitionedHypergraph(const PartitionID k, Hypergraph& hypergraph) :
    _k(k),
    _node(hypergraph.numaNode()),
    _hg(&hypergraph),
    part_weight(k, CAtomic<HypernodeWeight>(0)),
    part(hypergraph.initialNumNodes(), kInvalidPartition),
    pins_in_part(hypergraph.initialNumEdges() * k, CAtomic<HypernodeID>(0)),
    connectivity_sets(hypergraph.initialNumEdges(), k),
    move_to_penalty(hypergraph.initialNumNodes() * k, CAtomic<HyperedgeWeight>(0)),
    move_from_benefit(hypergraph.initialNumNodes() * k, CAtomic<HyperedgeWeight>(0))
  {

  }

  // TODO bring back the parallel memory allocation once the redesign is finalized
  explicit PartitionedHypergraph(const PartitionID k, const TaskGroupID, Hypergraph& hypergraph) : PartitionedHypergraph(k, hypergraph) { }

  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other) = default;
  PartitionedHypergraph & operator= (PartitionedHypergraph&& other) = default;

  ~PartitionedHypergraph() {
    freeInternalData();
  }

  bool setOnlyNodePart(const HypernodeID u, PartitionID p) {
    partAssertions(p);
    const HypernodeID u_local = common::get_local_position_of_vertex(u);
    ASSERT(part[u_local] == kInvalidPartition);
    part[u_local] = p;
    return true;
  }

  bool setNodePart(const HypernodeID u, PartitionID p) {
    if (setOnlyNodePart(u, p)) {
      part_weight[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
      for (HyperedgeID he : incidentEdges(u)) {
        incrementPinCountInPartWithoutGainUpdate(he, p);
      }
      return true;
    }
    return false;
  }

  bool setNodePart(const HypernodeID u, PartitionID p, parallel::scalable_vector<PartitionedHypergraph>& hypergraphs) {
    if (setOnlyNodePart(u, p)) {
      part_weight[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
      for (HyperedgeID he : incidentEdges(u)) {
        common::hypergraph_of_edge(he, hypergraphs).incrementPinCountInPartWithoutGainUpdate(he, p);
      }
      return true;
    }
    return false;
  }

  // This function type does not update part weights but instead updates a slack
  // If this is used, part weights have to be recomputed at the end of partition initialization
  // This function serves only as the counterpart to the more useful changeNodePartWithBalanceCheck function
  // The use of thread-local slack variables is supposed to alleviate contention on the part_weight vector
  bool setOnlyNodePartWithBalanceCheck(const HypernodeID u, PartitionID p, CAtomic<HypernodeWeight>& budget_p) {
    const HypernodeWeight wu = nodeWeight(u);
    if (budget_p.sub_fetch(wu, std::memory_order_relaxed) >= 0) {
      part[common::get_local_position_of_vertex(u)] = p;
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
    part_weight[to].fetch_add(wu, std::memory_order_relaxed);
    part_weight[from].fetch_sub(wu, std::memory_order_relaxed);
    part[common::get_local_position_of_vertex(u)] = to;
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  bool changeNodePart(const HypernodeID u, PartitionID from, PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    changeOnlyNodePart(u, from,  to);
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      HypernodeID pin_count_after_to = incrementPinCountInPartWithoutGainUpdate(he, to);
      HypernodeID pin_count_after_from = decrementPinCountInPartWithoutGainUpdate(he, from);
      delta_func(he, edgeWeight(he), edgeSize(he), pin_count_after_from, pin_count_after_to);
    }
    return true;
  }

  bool changeNodePart(const HypernodeID u, PartitionID from, PartitionID to,
                      parallel::scalable_vector<PartitionedHypergraph>& hgs,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
   changeOnlyNodePart(u, from, to);
   for ( const HyperedgeID& he : incidentEdges(u) ) {
      auto& hg_of_he = common::hypergraph_of_edge(he, hgs);
      HypernodeID pin_count_after_to = hg_of_he.incrementPinCountInPartWithoutGainUpdate(he, to);
      HypernodeID pin_count_after_from = hg_of_he.decrementPinCountInPartWithoutGainUpdate(he, from);
      delta_func(he, hg_of_he.edgeWeight(he), hg_of_he.edgeSize(he), pin_count_after_from, pin_count_after_to);
    }
    return true;
  }

  // Additionally rejects the requested move if it violates balance
  // must recompute part_weights once finished moving nodes
  bool changeNodePartWithBalanceCheckAndGainUpdates(const HypernodeID u,
                                                    PartitionID from,
                                                    CAtomic<HypernodeWeight>& budget_from,
                                                    PartitionID to,
                                                    CAtomic<HypernodeWeight>& budget_to) {
    const bool success = setOnlyNodePartWithBalanceCheck(u, to, budget_to);
    if (success) {
      budget_from.fetch_add(nodeWeight(u), std::memory_order_relaxed);
      for (HyperedgeID he : incidentEdges(u)) {
        incrementPinCountInPartWithGainUpdate(he, to);
        decrementPinCountInPartWithGainUpdate(he, from);
      }
    }
    return success;
  }

  void applyPartWeightUpdates(vec<HypernodeWeight>& part_weight_deltas) {
    for (PartitionID p = 0; p < _k; ++p) {
      part_weight[p].fetch_add(part_weight_deltas[p], std::memory_order_relaxed);
    }
  }

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    nodeAssertions(u);
    return part[common::get_local_position_of_vertex(u)];
  }

  void initializeBlockWeights() {
    auto accumulate = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HypernodeWeight> pws(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u_local = r.begin(); u_local < r.end(); ++u_local) {
        const PartitionID pu = partID( common::get_global_vertex_id(_node, u_local) );
        const HypernodeWeight wu = nodeWeight( common::get_global_vertex_id(_node, u_local) );
        pws[pu] += wu;
      }
      applyPartWeightUpdates(pws);
    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  template<typename FGetPartID>
  void initializePinCountInPart(FGetPartID get_part_id) {
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);

    auto assign = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HypernodeID>& pin_counts = ets_pin_count_in_part.local();
      for (HyperedgeID he_local = r.begin(); he_local < r.end(); ++he_local) {
        const HyperedgeID he_global = common::get_global_edge_id(_node, he_local);

        for (HypernodeID pin : pins(he_global)) {
          ++pin_counts[get_part_id(pin)];
        }

        for (PartitionID p = 0; p < _k; ++p) {
          assert(pinCountInPart(he_global, p) == 0);
          if (pin_counts[p] > 0) {
            connectivity_sets.add(he_local, p);
            pins_in_part[he_local * _k + p].store(pin_counts[p], std::memory_order_relaxed);
          }
          pin_counts[p] = 0;
        }

      }
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HyperedgeID(0), initialNumEdges()), assign);
  }


  // requires pinCountInPart to be computed
  void initializeGainInformation() {
    // check whether part has been initialized
    ASSERT( std::all_of(part.begin(), part.end(), [&](const PartitionID p) { return p != kInvalidPartition; }) );

    // we can either assume that pinCountInPart has been initialized and use this information
    // or recompute the gain information from scratch, which would be independent but requires more memory volume
    // since the LP refiner currently does not use and update the gain information, we rely on pinCountInPart being up to date
    // as we have to call initializeGainInformation() before every FM call
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_mfb(_k, 0), ets_mtp(_k, 0);

    auto accumulate_and_assign = [&](tbb::blocked_range<HypernodeID>& r) {

      vec<HyperedgeWeight>& l_move_from_benefit = ets_mfb.local();
      vec<HyperedgeWeight>& l_move_to_penalty = ets_mtp.local();

      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        for (HyperedgeID he : incidentEdges(u)) {
          HyperedgeWeight we = edgeWeight(he);
          for (PartitionID p = 0; p < _k; ++p) {
            const HypernodeID pcip = pinCountInPart(he, p);
            if (pcip == 0) {
              l_move_to_penalty[p] += we;
            } else if (pcip == 1) {
              l_move_from_benefit[p] += we;
            }
          }
        }

        for (PartitionID p = 0; p < _k; ++p) {
          move_to_penalty[u * _k + p].store(l_move_to_penalty[p], std::memory_order_relaxed);
          move_from_benefit[u * _k + p].store(l_move_from_benefit[p], std::memory_order_relaxed);
          // reset for next round
          l_move_to_penalty[p] = 0;
          l_move_from_benefit[p] = 0;
        }
      }

    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()), accumulate_and_assign);
  }

  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, part info, pin counts in part and border
  // ! vertices have to be computed in a postprocessing step.
  void initializePartition(const TaskGroupID ) {
    tbb::parallel_invoke(
            [&] { initializeBlockWeights(); },
            [&] { initializePinCountInPart( [&](HypernodeID u) { return partID(u); } ); }
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
    HypernodeID result = 0;
    for (const HyperedgeID he : incidentEdges(u)) {
      result += connectivity(he) > 1 ? 1 : 0;
    }
    return result;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    edgeAssertions(e);
    return connectivity_sets.connectivity(common::get_local_position_of_edge(e));
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    pinCountInPartAssertions(e, p);
    return pins_in_part[common::get_local_position_of_edge(e) * _k + p].load(std::memory_order_relaxed);
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID p) const {
    partAssertions(p);
    return part_weight[p];
  }

  HyperedgeWeight moveFromBenefit(const HypernodeID u, PartitionID p) const {
    return move_from_benefit[common::get_local_position_of_vertex(u) * _k + p].load(std::memory_order_relaxed);
  }

  HyperedgeWeight moveToPenalty(const HypernodeID u, PartitionID p) const {
    return move_to_penalty[common::get_local_position_of_vertex(u) * _k + p].load(std::memory_order_relaxed);
  }

  HyperedgeWeight km1Gain(const HypernodeID u, PartitionID from, PartitionID to) const {
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveFromBenefit(u, from) - moveToPenalty(u, to);
  }

  // TODO extend to a version with balance checks
  std::pair<PartitionID, HyperedgeWeight> bestDestinationBlock(HypernodeID u) const {
    HyperedgeWeight least_penalty = std::numeric_limits<HyperedgeWeight>::max();
    PartitionID best_index = 0;
    const PartitionID first = common::get_local_position_of_vertex(u) * _k;
    const PartitionID firstInvalid = first + _k;
    for (PartitionID index = first; index < firstInvalid; ++index) {
      const HyperedgeWeight penalty = move_to_penalty[index].load(std::memory_order_relaxed);
      if (penalty < least_penalty) {
        least_penalty = penalty;
        best_index = index;
      }
    }
    return std::make_pair(best_index - first, least_penalty);
  };

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    part.assign(part.size(), kInvalidPartition);
    for (auto& x : pins_in_part) x.store(0, std::memory_order_relaxed);
    for (auto& x : part_weight) x.store(0, std::memory_order_relaxed);
    connectivity_sets.reset();
  }

  auto& getPinCountInPartVector() {
    return pins_in_part;
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      tbb::parallel_invoke( [&] { parallel::free(pins_in_part); }, [&] { connectivity_sets.freeInternalData(); } );
    }
    _k = 0;
  }

 private:

  void edgeAssertions(const HyperedgeID e) const {
    unused(e);
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    const HyperedgeID local_id = common::get_local_position_of_edge(e);
    unused(local_id);
    ASSERT(local_id < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(_node == common::get_numa_node_of_edge(e), "Hyperedge" << e << "is not part of numa node" << _node);
  }

  void nodeAssertions(const HypernodeID u) const {
    ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled");
    const HyperedgeID local_id = common::get_local_position_of_vertex(u);
    unused(local_id);
    ASSERT(local_id < initialNumNodes(), "Hypernode" << u << "does not exist");
    if (_node != common::get_numa_node_of_vertex(u)) {
      LOG << V(u) << V(_node);
    }
    ASSERT(_node == common::get_numa_node_of_vertex(u), "Hypernode" << u << "is not part of numa node" << _node);
  }

  void partAssertions(const PartitionID p) const {
    unused(p);
    ASSERT(p != kInvalidPartition && p < _k);
  }

  void pinCountInPartAssertions(const HyperedgeID e, const PartitionID p) const {
    edgeAssertions(e);
    partAssertions(p);
  }

  void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
    nodeAssertions(u);
    partAssertions(p);
    ASSERT(common::get_local_position_of_vertex(u) * _k + p < move_from_benefit.size());
  }

  HypernodeID decrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    pinCountInPartAssertions(e, p);
    const HyperedgeID e_local = common::get_local_position_of_edge(e);
    const HypernodeID pin_count_after = pins_in_part[e_local * _k + p].sub_fetch(1, std::memory_order_relaxed);
    if ( pin_count_after == 0 ) {
      connectivity_sets.remove(e_local, p);
    }
    return pin_count_after;
  }

  HypernodeID incrementPinCountInPartWithoutGainUpdate(const HyperedgeID e, const PartitionID p) {
    pinCountInPartAssertions(e, p);
    const HyperedgeID e_local = common::get_local_position_of_edge(e);
    const HypernodeID pin_count_after = pins_in_part[e_local * _k + p].add_fetch(1, std::memory_order_relaxed);
    if ( pin_count_after == 1 ) {
      connectivity_sets.add(e_local, p);
    }
    return pin_count_after;
  }

  HypernodeID decrementPinCountInPartWithGainUpdate(const HyperedgeID e, const PartitionID p) {
    const HypernodeID pin_count_after = decrementPinCountInPartWithoutGainUpdate(e, p);
    if (pin_count_after == 1) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        const HypernodeID u_local = common::get_local_position_of_vertex(u);
        move_from_benefit[u_local * _k + p].fetch_add(we, std::memory_order_relaxed);
      }
    } else if (pin_count_after == 0) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        const HypernodeID u_local = common::get_local_position_of_vertex(u);
        move_from_benefit[u_local * _k + p].fetch_sub(we, std::memory_order_relaxed);
        move_to_penalty[u_local *_k + p].fetch_add(we, std::memory_order_relaxed);
      }
    }
    return pin_count_after;
  }

  HypernodeID incrementPinCountInPartWithGainUpdate(const HyperedgeID e, const PartitionID p) {
    const HypernodeID pin_count_after = incrementPinCountInPartWithoutGainUpdate(e, p);
    if (pin_count_after == 1) {
      const HyperedgeWeight we = edgeWeight(e);
      for (HypernodeID u : pins(e)) {
        nodeGainAssertions(u, p);
        const HypernodeID u_local = common::get_local_position_of_vertex(u);
        move_from_benefit[u_local * _k + p].fetch_add(we, std::memory_order_relaxed);
        move_to_penalty[u_local * _k + p].fetch_sub(we, std::memory_order_relaxed);
      }
    }
    return pin_count_after;
  }

  // ! Number of blocks
  PartitionID _k = 0;

  // ! NUMA node of this partitioned hypergraph
  int _node = 0;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Weight and information for all blocks.
  vec< CAtomic<HypernodeWeight> > part_weight;

  // ! Current block IDs of the vertices
  vec< PartitionID > part;

  // ! For each hyperedge and each block, _pins_in_part stores the number of pins in that block
  vec< CAtomic<HypernodeID> > pins_in_part;

  // TODO we probably don't need connectivity sets any more, except in the IP hypergraphs which don't need parallelism support
  // ! For each hyperedge, _connectivity_sets stores the set of blocks that the hyperedge spans
  ConnectivitySets connectivity_sets;

  // ! For each node and block, the sum of incident edge weights with zero pins in that part
  vec< CAtomic<HyperedgeWeight> > move_to_penalty;

  // ! For each node and block, the sum of incident edge weights with exactly one pin in that part
  vec< CAtomic<HyperedgeWeight> > move_from_benefit;

public:

  // ! Returns a range to loop over the set of block ids contained in hyperedge e.
  IteratorRange<ConnectivitySets::Iterator> connectivitySet(const HyperedgeID e) const {
    edgeAssertions(e);
    return connectivity_sets.connectivitySet(common::get_local_position_of_edge(e));
  }


  // ! Restores a hyperedge of a certain size.
  void restoreEdge(const HyperedgeID he, const size_t size,
                   const HyperedgeID representative = kInvalidHyperedge) {
    _hg->restoreEdge(he, size, representative);
    for ( const HypernodeID& pin : pins(he) ) {
      incrementPinCountInPartWithoutGainUpdate(he, partID(pin));
      // this doesn't break km1 gains, but unfortunately it does break the actual penalty and benefit values
      // we would have to add edgeWeight(he) to all move_to_penalty[pin * _k + p] for p != partID(pin)
      // and add it to move_from_benefit[pin * _k + partID(pin)]
      // I'm assuming this is only used in postprocessing to restore removed single pin hyperedges in the preprocessing?
      // Then we don't need those values.
      // Apply same decision in numa_partitioned_hypergraph
      // TODO decide which version we want
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


  // Hypergraph Forwards

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

  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
// TODO finish this function when everything else is done

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);
    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    connectivity_sets.memoryConsumption(connectivity_set_node);
    parent->addChild("Pin Count In Part", sizeof(CAtomic<HypernodeID>) * _k * _hg->initialNumEdges());
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


};

} // namespace ds
} // namespace mt_kahypar
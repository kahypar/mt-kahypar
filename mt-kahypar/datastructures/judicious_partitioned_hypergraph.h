/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2022 Noah Wahl <noah.wahl@student.kit.edu>
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
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar {
namespace ds {

class JudiciousPartitionedHypergraph {
private:
  using Hypergraph = StaticHypergraph;
  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

public:
  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;

  JudiciousPartitionedHypergraph() = default;

  explicit JudiciousPartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
    _k(k),
    _hg(&hypergraph),
    _part_loads(k, 0),
    _part_ids(hypergraph.initialNumNodes(), kInvalidPartition),
    _pins_in_part(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize(), false),
    _connectivity_set(hypergraph.initialNumEdges(), k, false) {
  }

  JudiciousPartitionedHypergraph(const JudiciousPartitionedHypergraph&) = delete;
  JudiciousPartitionedHypergraph & operator= (const JudiciousPartitionedHypergraph &) = delete;

  JudiciousPartitionedHypergraph(JudiciousPartitionedHypergraph&& other) = default;
  JudiciousPartitionedHypergraph & operator= (JudiciousPartitionedHypergraph&& other) = default;

  ~JudiciousPartitionedHypergraph() {
    freeInternalData();
  }

  void resetData() {
    tbb::parallel_invoke([&] {
      _part_ids.assign(_part_ids.size(), kInvalidPartition);
    }, [&] {
      _pins_in_part.data().assign(_pins_in_part.data().size(), 0);
    }, [&] {
      _connectivity_set.reset();
    }, [&] {
      _part_loads.assign(_part_loads.size(), 0);
    });
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
    static_cast<const JudiciousPartitionedHypergraph&>(*this).doParallelForAllNodes(f);
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
    static_cast<const JudiciousPartitionedHypergraph&>(*this).doParallelForAllEdges(f);
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

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u, const PartitionID to) {
    _hg->restoreDegreeZeroHypernode(u);
    setNodePart(u, to);
  }

  HyperedgeWeight weightOfDisabledEdges(const HypernodeID n) const {
    return _hg->weightOfDisabledEdges(n);
  }


  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return _hg->edgeWeight(e);
  }

  // ! Unique id of a hyperedge
  HyperedgeID uniqueEdgeID(const HyperedgeID e) const {
    return e;
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

  // ####################### Partition Information #######################

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    return _part_ids[u];
  }

  void extractPartIDs(vec<PartitionID>& part_ids) {
    std::swap(_part_ids, part_ids);
  }

  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(u < _part_ids.size());
    ASSERT(_part_ids[u] == kInvalidPartition);
    _part_ids[u] = p;
  }

  void setNodePart(const HypernodeID u, PartitionID p) {
    setOnlyNodePart(u, p);
    ASSERT(p < _k);
    _part_loads[p] += weightOfDisabledEdges(u);
    for (HyperedgeID he : incidentEdges(u)) {
      ASSERT(edgeIsEnabled(he));
      const HypernodeID pin_count_after = incrementPinCountInPartWithoutGainUpdate(he, p);
      if (pin_count_after == 1) {
        _part_loads[p] += edgeWeight(he);
      }
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.

  template<typename SuccessFunc, typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      SuccessFunc&& report_success,
                      DeltaFunc&& delta_func) {
    assert(partID(u) == from);
    assert(from != to);
    ASSERT(u < _part_ids.size());
    _part_ids[u] = to;
    ASSERT(_part_loads[from] > 0);
    _part_loads[from] -= weightOfDisabledEdges(u);
    _part_loads[to] += weightOfDisabledEdges(u);
    report_success();
    for ( const HyperedgeID he : incidentEdges(u) ) {
      ASSERT(edgeIsEnabled(he));
      updatePinCountOfHyperedge(he, from, to, delta_func);
    }
    return true;
  }

  // curry
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC) {
    return changeNodePart(u, from, to, []{}, delta_func);
  }

  // ! Load of a block
  HyperedgeWeight partLoad(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_loads[p];
  }

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    HypernodeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( connectivity(he) > 1 ) {
        ++num_incident_cut_hyperedges;
      }
    }
    return num_incident_cut_hyperedges;
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

  void initializePartition() {
    initializePinCountInPart();
    initializeBlockLoads();
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    _part_ids.assign(_part_ids.size(), kInvalidPartition);
    _part_loads.assign(_part_loads.size(), kInvalidPartition);

    // Reset pin count in part and connectivity set
    for ( const HyperedgeID& he : edges() ) {
      for ( const PartitionID& block : connectivitySet(he) ) {
        _pins_in_part.setPinCountInPart(he, block, 0);
      }
      _connectivity_set.clear(he);
    }
  }


  // ####################### Memory Consumption #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);
    utils::MemoryTreeNode* connectivity_set_node = parent->addChild("Connectivity Sets");
    _connectivity_set.memoryConsumption(connectivity_set_node);

    parent->addChild("Part Loads", sizeof(HyperedgeWeight) * _k);
    parent->addChild("Part IDs", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("Pin Count In Part", _pins_in_part.size_in_bytes());
  }

  void freeInternalData() {
    if ( _k > 0 ) {
      _connectivity_set.freeInternalData();
    }
    _k = 0;
  }

private:

  void applyPartLoadUpdates(vec<HyperedgeWeight>& part_load_deltas) {
    for (PartitionID p = 0; p < _k; ++p) {
      __atomic_fetch_add(&_part_loads[p], part_load_deltas[p], __ATOMIC_SEQ_CST);
    }
  }

  void initializeBlockLoads() {
    auto accumulate = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HyperedgeWeight> pvs(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HyperedgeID e = r.begin(); e < r.end(); ++e) {
        if ( edgeIsEnabled(e) ) {
          for (const auto p : connectivitySet(e)) {
            pvs[p] += edgeWeight(e);
          }
        }
      }
      applyPartLoadUpdates(pvs);
    };
    auto accumulate_disabled_edges = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HyperedgeWeight> pvs(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if (nodeIsEnabled(u)) {
          pvs[partID(u)] += weightOfDisabledEdges(u);
        }
      }
      applyPartLoadUpdates(pvs);
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HypernodeID(0), initialNumEdges()),
                      accumulate,
                      tbb::static_partitioner()
    );
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate_disabled_edges,
                      tbb::static_partitioner()
    );
  }

  void initializePinCountInPart() {
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);

    auto assign = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HypernodeID>& pin_counts = ets_pin_count_in_part.local();
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {
        if ( edgeIsEnabled(he) ) {
          for (const HypernodeID& pin : pins(he)) {
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

  // ! Updates pin count in part using a spinlock.
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePinCountOfHyperedge(const HyperedgeID he,
                                                                    const PartitionID from,
                                                                    const PartitionID to,
                                                                    const DeltaFunction& delta_func) {
    const HypernodeID pin_count_in_from_part_after = decrementPinCountInPartWithoutGainUpdate(he, from);
    const HypernodeID pin_count_in_to_part_after = incrementPinCountInPartWithoutGainUpdate(he, to);
    if (pin_count_in_from_part_after == 0) {
      ASSERT(_part_loads[from] > 0);
      _part_loads[from] -= edgeWeight(he);
    }
    if (pin_count_in_to_part_after == 1) {
      _part_loads[to] += edgeWeight(he);
    }
    delta_func(he, edgeWeight(he), edgeSize(he), pin_count_in_from_part_after, pin_count_in_to_part_after);
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


private:
  // ! Number of blocks
  PartitionID _k = 0;

  // ! Hypergraph object around which this partitioned hypergraph is wrapped
  Hypergraph* _hg = nullptr;

  // ! Load and information for all blocks.
  vec< HyperedgeWeight > _part_loads;

  // ! Current block IDs of the vertices
  vec< PartitionID > _part_ids;

  // ! For each hyperedge and each block, _pins_in_part stores the
  // ! number of pins in that block
  PinCountInPart _pins_in_part;

  // ! For each hyperedge, _connectivity_set stores the set of blocks that the hyperedge spans
  ConnectivitySets _connectivity_set;
};

}
}

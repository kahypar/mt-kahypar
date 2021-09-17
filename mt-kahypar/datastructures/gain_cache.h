//
// Created by mlaupichler on 02.06.21.
//

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/async/iterable_bit_set.h"
#include "mt-kahypar/datastructures/async/no_downsize_integral_type_vector.h"

namespace mt_kahypar::ds {

    using CompressedConnectivitySetSnapshot = NoDownsizeIntegralTypeVector<PartitionID>;

    using PartIDs = Array<PartitionID>;

class LightGainCache {

public:

    LightGainCache() = default;

    LightGainCache(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _part_ids(partIDs),
        _pins_in_parts(pin_count_in_part),
        _connectivity_sets(connectivity_sets),
        _move_from_benefit(),
        _move_to_penalty() {}

    LightGainCache(const HypernodeID num_nodes, const PartitionID k,
                     const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) :
        _num_nodes(num_nodes),
        _k(k),
        _part_ids(partIDs),
        _pins_in_parts(pin_count_in_part),
        _connectivity_sets(connectivity_sets),
        _move_from_benefit("Refinement", "benefit_cache", num_nodes, true, false),
        _move_to_penalty("Refinement", "penalty_cache", size_t(num_nodes) * size_t(k + 1), true, false) {}


    void assignQueryObjects(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) {
      _part_ids = partIDs;
      _pins_in_parts = pin_count_in_part;
      _connectivity_sets = connectivity_sets;
    }

    void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) {
      ASSERT(_num_nodes == 0 && _k == kInvalidPartition);
      _num_nodes = num_nodes;
      _k = k;
      _move_from_benefit.resize(group, key, num_nodes, true, true);
      _move_to_penalty.resize(group, key, size_t(num_nodes) * size_t(k + 1), true, true);
    }

    size_t size_in_bytes() const {
      return sizeof(CAtomic<HyperedgeWeight>) * (_move_from_benefit.size() + _move_to_penalty.size());
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
      return _move_from_benefit[u].load(m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID, const PartitionID, std::memory_order m = std::memory_order_relaxed) const {
      unused(m);
      ERROR("moveFromBenefit(HypernodeID, PartitionID) not supported on LightGainCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
      return _move_to_penalty[incident_net_weight_index(u)].load(m) -
             _move_to_penalty[penalty_index(u, to)].load(m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator, std::memory_order m = std::memory_order_relaxed) {
      ASSERT(_part_ids);
      auto block_of_u = (*_part_ids)[u];
      _move_from_benefit[u].store(benefit_aggregator[block_of_u], m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator,
                         const HyperedgeWeight incident_edges_weight, vec<HyperedgeWeight>& penalty_aggregator) {
      ASSERT(_part_ids);

      // Benefit
      auto block_of_u = (*_part_ids)[u];
      _move_from_benefit[u].store(benefit_aggregator[block_of_u], std::memory_order_relaxed);
      // Reset all entries to zero
      benefit_aggregator.assign(_k,0);

      // Penalty
      _move_to_penalty[incident_net_weight_index(u)].store(
          incident_edges_weight, std::memory_order_relaxed);
      for (PartitionID p = 0; p < _k; ++p) {
        _move_to_penalty[penalty_index(u,p)].store(
            penalty_aggregator[p], std::memory_order_relaxed);
        penalty_aggregator[p] = 0;
      }
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v,
                                        const PartitionID block,
                                        const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins) {
      ASSERT(_part_ids);

      if ( pin_count_in_part_after == 2 ) {
        // u might be replaced by an other vertex in the batch
        // => search for other pin of the corresponding block and
        // subtract edge weight.
        for ( auto it = pins.begin(); it != pins.end(); ++it ) {
          const HypernodeID& pin = *it;
          if ( pin != v && (*_part_ids)[pin] == block ) {
            _move_from_benefit[pin].sub_fetch(we, std::memory_order_relaxed);
            break;
          }
        }
      }

      ASSERT(_connectivity_sets);
      // For all blocks contained in the connectivity set of hyperedge he
      // we increase the move_to_penalty for vertex v by w(e) =>
      // move_to_penalty is then w(I(v)) - move_to_penalty(v, p) for a
      // block p
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const PartitionID block_in_conn_set : _connectivity_sets->connectivitySet(he)) {
        _move_to_penalty[penalty_index(v, block_in_conn_set)].add_fetch(
            we, std::memory_order_relaxed);
      }
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight, const HypernodeID, const PartitionID,
                                         const HypernodeID, IteratorRange<PinIteratorT>,
                                         const ConnectivitySets::Snapshot&, const PinCountInPart::Snapshot&) {
      ERROR("asyncUpdateForUncontractCaseOne(HyperedgeWeight, HypernodeID, PartitionID, HypernodeID, "
            "IteratorRange<PinIteratorT>, const ConnectivitySets::Snapshot&, const PinCountInPart::Snapshot&) not supported on LightGainCache.");
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight, const HypernodeID, const PartitionID,
                                         const HypernodeID, IteratorRange<PinIteratorT>,
                                         const CompressedConnectivitySetSnapshot &, const CompressedConnectivitySetSnapshot &) {
      ERROR("asyncUpdateForUncontractCaseOne(HyperedgeWeight, HypernodeID, PartitionID, HypernodeID, "
            "IteratorRange<PinIteratorT>, const CompressedConnectivitySetSnapshot &, const CompressedConnectivitySetSnapshot &) not supported on LightGainCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isPinCountThatTriggersUpdateForAllPinsForUncontractCaseOne(const HypernodeID) {
      ERROR("isPinCountThatTriggersUpdateForAllPinsForUncontractCaseOne() should not be called on LightGainCache.");
    }


    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u, const HypernodeID v) {
      // Since u is no longer incident to hyperedge he its contribution for decreasing
      // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).
      ASSERT(_part_ids);
      ASSERT(_pins_in_parts);
      auto block_of_u = (*_part_ids)[u];
      if ( _pins_in_parts->pinCountInPart(he, block_of_u) == 1 ) {
        _move_from_benefit[u].sub_fetch(we, std::memory_order_relaxed);
        _move_from_benefit[v].add_fetch(we, std::memory_order_relaxed);
      }

      ASSERT(_connectivity_sets);
      // For all blocks contained in the connectivity set of hyperedge he
      // we decrease the the move_to_penalty for vertex u and increase it for
      // vertex v by w(e)
      _move_to_penalty[incident_net_weight_index(u)].sub_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const PartitionID block : _connectivity_sets->connectivitySet(he)) {
        _move_to_penalty[penalty_index(u, block)].sub_fetch(
            we, std::memory_order_relaxed);
        _move_to_penalty[penalty_index(v, block)].add_fetch(
            we, std::memory_order_relaxed);
      }
    }

        MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight, const HypernodeID, const HypernodeID,
                                         const ConnectivitySets::Snapshot&, const PinCountInPart::Snapshot&) {
      ERROR("asyncUpdateForUncontractCaseTwo(HyperedgeWeight, HypernodeID, HypernodeID, HypernodeID,"
            "const ConnectivitySets::Snapshot&, const PinCountInPart::Snapshot&) not supported on LightGainCache.");
    }

    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight, const HypernodeID, const HypernodeID,
                                         const CompressedConnectivitySetSnapshot &, const CompressedConnectivitySetSnapshot &) {
      ERROR("asyncUpdateForUncontractCaseTwo(HyperedgeWeight, HypernodeID, HypernodeID, HypernodeID,"
            "const CompressedConnectivitySetSnapshot &, const CompressedConnectivitySetSnapshot &) not supported on LightGainCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block_of_single_pin) {
      _move_from_benefit[single_pin_in_he].add_fetch(
          we, std::memory_order_relaxed);

      _move_to_penalty[incident_net_weight_index(single_pin_in_he)].add_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[penalty_index(single_pin_in_he, block_of_single_pin)].add_fetch(
          we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
      ASSERT(_part_ids);
      if (pin_count_in_from_part_after == 0) {
        for (const auto& u : pins) {
          nodeGainAssertions(u, from);
          _move_to_penalty[penalty_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
        }
      } else if (pin_count_in_from_part_after == 1) {
        for (const auto& u :pins) {
          nodeGainAssertions(u, from);
          if ((*_part_ids)[u] == from) {
            _move_from_benefit[u].fetch_add(we, std::memory_order_relaxed);
          }
        }
      }

      if (pin_count_in_to_part_after == 1) {
        for (const auto& u : pins) {
          nodeGainAssertions(u, to);
          _move_to_penalty[penalty_index(u, to)].fetch_add(we, std::memory_order_relaxed);
        }
      } else if (pin_count_in_to_part_after == 2) {
        for (const auto& u :pins) {
          nodeGainAssertions(u, to);
          if ((*_part_ids)[u] == to) {
            _move_from_benefit[u].fetch_sub(we, std::memory_order_relaxed);
          }
        }
      }
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool arePinCountsAfterMoveThatTriggerUpdateForAllPins(const HypernodeID, const HypernodeID) {
      ERROR("arePinCountsAfterMoveThatTriggerUpdateForAllPins should not be called on LightGainCache.");
    }

private:

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidPenaltyEntry(const HypernodeID u, const PartitionID p) const {
      return penalty_index(u, p) < _move_to_penalty.size();
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    size_t incident_net_weight_index(const HypernodeID u) const {
      return size_t(u) * ( _k + 1 );
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    size_t penalty_index(const HypernodeID u, const PartitionID p) const {
      return size_t(u) * ( _k + 1 )  + p + 1;
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
      unused(u);
      unused(p);
      ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
//        ASSERT(_hg_query_funcs->is_node_enabled(u), "Hypernode" << u << "is disabled");
      ASSERT(p != kInvalidPartition && p < _k);
      ASSERT(isValidPenaltyEntry(u,p));
    }

    // ! Total number of hypernodes that gains are stored for
    HypernodeID _num_nodes;

    // ! Number of partitions that penalties are stored for per hypernode
    PartitionID _k;

    const PartIDs* _part_ids;
    const PinCountInPart* _pins_in_parts;
    const ConnectivitySets* _connectivity_sets;

    Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

    // ! For each node and block, the sum of weights of incident edges with zero pins in that block
    Array< CAtomic<HyperedgeWeight> > _move_to_penalty;
};

class HeavyGainCache {

    public:

    HeavyGainCache() = default;

    HeavyGainCache(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets, parallel_tag_t) :
            _num_nodes(0),
            _k(kInvalidPartition),
            _part_ids(partIDs),
            _pins_in_parts(pin_count_in_part),
            _connectivity_sets(connectivity_sets),
            _move_from_benefit(),
            _move_to_penalty() {}

    HeavyGainCache(const HypernodeID num_nodes, const PartitionID k,
                       const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) :
            _num_nodes(num_nodes),
            _k(k),
            _part_ids(partIDs),
            _pins_in_parts(pin_count_in_part),
            _connectivity_sets(connectivity_sets),
            _move_from_benefit("Refinement", "benefit_cache", static_cast<size_t>(num_nodes) * static_cast<size_t>(k), true, false),
            _move_to_penalty("Refinement", "penalty_cache", size_t(num_nodes) * size_t(k + 1), true, false) {}


        void assignQueryObjects(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) {
          _part_ids = partIDs;
          _pins_in_parts = pin_count_in_part;
          _connectivity_sets = connectivity_sets;
        }

        void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) {
          ASSERT(_num_nodes == 0 && _k == kInvalidPartition);
          _num_nodes = num_nodes;
          _k = k;
          _move_from_benefit.resize(group, key, static_cast<size_t>(num_nodes) * static_cast<size_t>(k), true, true);
          _move_to_penalty.resize(group, key, size_t(num_nodes) * size_t(k + 1), true, true);
        }

        size_t size_in_bytes() const {
          return sizeof(CAtomic<HyperedgeWeight>) * (_move_from_benefit.size() + _move_to_penalty.size());
        }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
      ASSERT(_part_ids);
      return _move_from_benefit[benefit_index(u, (*_part_ids)[u])].load(m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, const PartitionID p, std::memory_order m = std::memory_order_relaxed) const {
      return _move_from_benefit[benefit_index(u, p)].load(m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator, std::memory_order m = std::memory_order_relaxed) {
      for (PartitionID i = 0; i < _k; ++i) {
        _move_from_benefit[benefit_index(u,i)].store(benefit_aggregator[i], m);
      }
    }

        MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
        HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
          return _move_to_penalty[incident_net_weight_index(u)].load(m) -
                 _move_to_penalty[penalty_index(u, to)].load(m);
        }

        MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
        void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator,
                             const HyperedgeWeight incident_edges_weight, vec<HyperedgeWeight>& penalty_aggregator) {

          _move_to_penalty[incident_net_weight_index(u)].store(
              incident_edges_weight, std::memory_order_relaxed);

          for (PartitionID p = 0; p < _k; ++p) {
            _move_from_benefit[benefit_index(u,p)].store(benefit_aggregator[p], std::memory_order_relaxed);
            _move_to_penalty[penalty_index(u,p)].store(penalty_aggregator[p], std::memory_order_relaxed);
            // Reset entry to zero so aggregator can be reused
            benefit_aggregator[p] = 0;
            penalty_aggregator[p] = 0;
          }
        }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> connectivity set queried on demand
    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v,
                                        const PartitionID block, const HypernodeID pin_count_in_part_after,
                                        IteratorRange<PinIteratorT> pins) {
      // Calculate and add contribution of he to the benefit of the uncontracted node v.
      // For all blocks contained in the connectivity set of hyperedge he
      // we increase the move_to_penalty for vertex v by w(e) =>
      // move_to_penalty is then w(I(v)) - move_to_penalty(v, p) for a
      // block p
      ASSERT(_connectivity_sets);
      ASSERT(_pins_in_parts);
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const auto & p : _connectivity_sets->connectivitySet(he)) {
        _move_to_penalty[penalty_index(v, p)].add_fetch(we, std::memory_order_relaxed);
        if (_pins_in_parts->pinCountInPart(he, p) == 1) {
          _move_from_benefit[benefit_index(v, p)].add_fetch(we, std::memory_order_relaxed);
        }
      }

      // In this case, u and v are incident to hyperedge he after uncontraction
      if ( pin_count_in_part_after == 2 ) {
        // Reduce benefits for other pins and this block (via snapshot of pins at that point), excluding v as it
        // never got the benefit to begin with
        for (const auto& pin : pins) {
          if (pin != v) {
            _move_from_benefit[benefit_index(pin, block)].sub_fetch(we, std::memory_order_relaxed);
          }
        }
      }
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v, const PartitionID block,
                                         const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins,
                                         const ConnectivitySets::Snapshot& connectivity_set_snapshot, const PinCountInPart::Snapshot& pin_count_in_part_snapshot) {
      // Calculate and add contribution of he to the benefit of the uncontracted node v
      // Add benefit for all blocks that include exactly one pin of he.
      // For all blocks contained in the connectivity set of hyperedge he
      // we increase the move_to_penalty for vertex v by w(e) =>
      // move_to_penalty is then w(I(v)) - move_to_penalty(v, p) for a
      // block p.
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const PartitionID block_in_conn_set : connectivity_set_snapshot.snapshottedConnectivitySet()) {
        _move_to_penalty[penalty_index(v, block_in_conn_set)].add_fetch(
            we, std::memory_order_relaxed);
        if (pin_count_in_part_snapshot.pinCountInPart(block_in_conn_set) == 1) {
          _move_from_benefit[benefit_index(v, block_in_conn_set)].add_fetch(we, std::memory_order_relaxed);
        }
      }

      // In this case, u and v are incident to hyperedge he after uncontraction
      if ( pin_count_in_part_after == 2 ) {
        // Reduce benefits for other pins and this block (via snapshot of pins at that point), excluding v as it
        // never got the benefit to begin with
        for (const auto& pin : pins) {
          if (pin != v) {
            _move_from_benefit[benefit_index(pin, block)].sub_fetch(we, std::memory_order_relaxed);
          }
        }
      }
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v, const PartitionID block,
                                         const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins,
                                         const CompressedConnectivitySetSnapshot & connectivity_set, const CompressedConnectivitySetSnapshot & parts_with_one_pin) {
      // Calculate and add contribution of he to the benefit of the uncontracted node v
      // Add benefit for all blocks that include exactly one pin of he.
      for (const auto& block_with_one_pin : parts_with_one_pin) {
        _move_from_benefit[benefit_index(v, block_with_one_pin)].add_fetch(we, std::memory_order_relaxed);
      }

      // In this case, u and v are incident to hyperedge he after uncontraction
      if ( pin_count_in_part_after == 2 ) {
        // Reduce benefits for other pins and this block (via snapshot of pins at that point), excluding v as it
        // never got the benefit to begin with
        for (const auto& pin : pins) {
          if (pin != v) {
            _move_from_benefit[benefit_index(pin, block)].sub_fetch(we, std::memory_order_relaxed);
          }
        }
      }

      // For all blocks contained in the connectivity set of hyperedge he
      // we increase the move_to_penalty for vertex v by w(e) =>
      // move_to_penalty is then w(I(v)) - move_to_penalty(v, p) for a
      // block p.
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const PartitionID block_in_conn_set : connectivity_set) {
        _move_to_penalty[penalty_index(v, block_in_conn_set)].add_fetch(
            we, std::memory_order_relaxed);
      }
    }

    // ! Allows user of gain cache to query whether the given pin count in a block after the uncontraction in case one
    // ! (u and v incident to edge after uncontraction) triggers an update of all pins of the edge. (Essentially returns whether the
    // ! pin count has been increased to two, demanding the PartitionedHypergraph to take a pin snapshot).
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isPinCountThatTriggersUpdateForAllPinsForUncontractCaseOne(const HypernodeID pin_count_in_part_after) {
      return pin_count_in_part_after == 2;
    }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> connectivity set queried on demand
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u,
                                        const HypernodeID v) {
      // Remove benefits contributed by this hyperedge from u as it no longer belongs to the hyperedge and add them to v.
      // For all blocks contained in the connectivity set of hyperedge he
      // we decrease the the move_to_penalty for vertex u and increase it for
      // vertex v by w(e)
      ASSERT(_connectivity_sets);
      ASSERT(_pins_in_parts);
      _move_to_penalty[incident_net_weight_index(u)].sub_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const auto & p : _connectivity_sets->connectivitySet(he)) {
        _move_to_penalty[penalty_index(u, p)].sub_fetch(
            we, std::memory_order_relaxed);
        _move_to_penalty[penalty_index(v, p)].add_fetch(
            we, std::memory_order_relaxed);
        if (_pins_in_parts->pinCountInPart(he, p) == 1) {
          _move_from_benefit[benefit_index(u, p)].sub_fetch(we, std::memory_order_relaxed);
          _move_from_benefit[benefit_index(v, p)].add_fetch(we, std::memory_order_relaxed);
        }
      }
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v,
                                         const ConnectivitySets::Snapshot& connectivity_set_snapshot, const PinCountInPart::Snapshot& pin_count_in_part_snapshot) {
      // In this case, u is replaced by v in hyperedge he
      // => Pin counts of hyperedge he does not change
      // Since u is no longer incident to hyperedge he its contribution for decreasing
      // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).

      // Remove benefits contributed by this hyperedge from u as it no longer belongs to the hyperedge and add them to v.
      // For all blocks contained in the connectivity set of hyperedge he
      // we decrease the the move_to_penalty for vertex u and increase it for
      // vertex v by w(e).
      _move_to_penalty[incident_net_weight_index(u)].sub_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const PartitionID block : connectivity_set_snapshot.snapshottedConnectivitySet()) {
        _move_to_penalty[penalty_index(u, block)].sub_fetch(
            we, std::memory_order_relaxed);
        _move_to_penalty[penalty_index(v, block)].add_fetch(
            we, std::memory_order_relaxed);
        if (pin_count_in_part_snapshot.pinCountInPart(block) == 1) {
          _move_from_benefit[benefit_index(u, block)].sub_fetch(we, std::memory_order_relaxed);
          _move_from_benefit[benefit_index(v, block)].add_fetch(we, std::memory_order_relaxed);
        }
      }
    }
    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v,
                                         const CompressedConnectivitySetSnapshot & connectivity_set, const CompressedConnectivitySetSnapshot & parts_with_one_pin) {
      // In this case, u is replaced by v in hyperedge he
      // => Pin counts of hyperedge he does not change
      // Since u is no longer incident to hyperedge he its contribution for decreasing
      // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).

      // Remove benefits contributed by this hyperedge from u as it no longer belongs to the hyperedge and add them to v.
      for (const auto& part_with_one_pin : parts_with_one_pin) {
        _move_from_benefit[benefit_index(u, part_with_one_pin)].sub_fetch(we, std::memory_order_relaxed);
        _move_from_benefit[benefit_index(v, part_with_one_pin)].add_fetch(we, std::memory_order_relaxed);
      }

      // For all blocks contained in the connectivity set of hyperedge he
      // we decrease the the move_to_penalty for vertex u and increase it for
      // vertex v by w(e).
      _move_to_penalty[incident_net_weight_index(u)].sub_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[incident_net_weight_index(v)].add_fetch(
          we, std::memory_order_relaxed);
      for (const PartitionID block : connectivity_set) {
        _move_to_penalty[penalty_index(u, block)].sub_fetch(
            we, std::memory_order_relaxed);
        _move_to_penalty[penalty_index(v, block)].add_fetch(
            we, std::memory_order_relaxed);
      }
    }


    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block_of_single_pin) {
      _move_from_benefit[benefit_index(single_pin_in_he, block_of_single_pin)].add_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[incident_net_weight_index(single_pin_in_he)].add_fetch(
          we, std::memory_order_relaxed);
      _move_to_penalty[penalty_index(single_pin_in_he, block_of_single_pin)].add_fetch(
          we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
      if (pin_count_in_from_part_after == 0) {
        for (const auto& u : pins) {
          nodeGainAssertions(u, from);
          _move_from_benefit[benefit_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
          _move_to_penalty[penalty_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
        }
      } else if (pin_count_in_from_part_after == 1) {
        for (const auto& u : pins) {
          nodeGainAssertions(u, from);
          _move_from_benefit[benefit_index(u, from)].fetch_add(we, std::memory_order_relaxed);
        }
      }

      if (pin_count_in_to_part_after == 1) {
        for (const auto& u : pins) {
          nodeGainAssertions(u, to);
          _move_from_benefit[benefit_index(u, to)].fetch_add(we, std::memory_order_relaxed);
          _move_to_penalty[penalty_index(u, to)].fetch_add(we, std::memory_order_relaxed);
        }
      } else if (pin_count_in_to_part_after == 2) {
        for (const auto& u : pins) {
          nodeGainAssertions(u, to);
          _move_from_benefit[benefit_index(u, to)].fetch_sub(we, std::memory_order_relaxed);
        }
      }
    }

    // ! Allows user of gain cache to query whether the given pin count in the origin or target block of an edge incident
    // ! to a moved node triggers an update of all pins of the edge. (Essentially returns whether the
    // ! pin count in the from part has been reduced to 0 or 1 or the pin count in the to part has been increased
    // ! to 1 or 2 demanding the PartitionedHypergraph to take a pin snapshot).
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool arePinCountsAfterMoveThatTriggerUpdateForAllPins(const HypernodeID pin_count_in_from_part_after, const HypernodeID pin_count_in_to_part_after) {
      return (pin_count_in_from_part_after == 0) || (pin_count_in_from_part_after == 1) || (pin_count_in_to_part_after == 1) || (pin_count_in_to_part_after == 2);
    }

    private:

        MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
        bool isValidEntry(const HypernodeID u, const PartitionID p) const {
          return penalty_index(u, p) < _move_to_penalty.size() && benefit_index(u, p) < _move_from_benefit.size();
        }

        MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
        size_t incident_net_weight_index(const HypernodeID u) const {
          return size_t(u) * ( _k + 1 );
        }

        MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
        size_t penalty_index(const HypernodeID u, const PartitionID p) const {
          return size_t(u) * ( _k + 1 )  + p + 1;
        }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    size_t benefit_index(const HypernodeID u, const PartitionID p) const {
      return size_t(u) * size_t(_k)  + size_t(p);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
      unused(u);
      unused(p);
      ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
//        ASSERT(_hg_query_funcs->is_node_enabled(u), "Hypernode" << u << "is disabled");
      ASSERT(p != kInvalidPartition && p < _k);
      ASSERT(isValidEntry(u,p));
    }

        // ! Total number of hypernodes that gains are stored for
        HypernodeID _num_nodes;

        // ! Number of partitions that penalties are stored for per hypernode
        PartitionID _k;

        const PartIDs* _part_ids;
        const PinCountInPart* _pins_in_parts;
        const ConnectivitySets* _connectivity_sets;

        // ! For each node and block, the sum of weights of incident edges with exactly one pin in that block
        Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

        // ! For each node and block, the sum of weights of incident edges with zero pins in that block
        Array< CAtomic<HyperedgeWeight> > _move_to_penalty;
    };

} // namespace mt_kahypar::ds
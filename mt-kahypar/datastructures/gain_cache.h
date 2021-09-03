//
// Created by mlaupichler on 02.06.21.
//

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/async/iterable_bit_set.h"
#include "mt-kahypar/datastructures/async/no_downsize_integral_type_vector.h"

namespace mt_kahypar::ds {

    using ConnectivitySetSnapshot = NoDownsizeIntegralTypeVector<PartitionID>;

    using PartIDs = Array<PartitionID>;

template<class MoveFromBenefitCache = Mandatory, class MoveToPenaltyCache = Mandatory>
class GainCacheFacade {

public:

    // Empty constructor to satisfy no-parameter constructor of PartitionedHypergraph
    GainCacheFacade() = default;

    // Default empty constructor to get size zero gain cache to allow subsequent parallel construction with parallel_resize()
    explicit GainCacheFacade(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _benefit_cache(partIDs, pin_count_in_part, connectivity_sets, parallel_tag_t()),
        _penalty_cache(partIDs, pin_count_in_part, connectivity_sets, parallel_tag_t()) {}

    GainCacheFacade(const HypernodeID num_nodes, const PartitionID k, const PartIDs* partIDs,
                    const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) :
        _num_nodes(num_nodes),
        _k(k),
        _benefit_cache(std::string(GROUP_NAME),std::string(BENEFIT_KEY), num_nodes, k, partIDs, pin_count_in_part, connectivity_sets),
        _penalty_cache(std::string(GROUP_NAME),std::string(BENEFIT_KEY), num_nodes, k, partIDs, pin_count_in_part, connectivity_sets) {}

    GainCacheFacade(const GainCacheFacade&) = delete;
    GainCacheFacade & operator= (const GainCacheFacade &) = delete;
    GainCacheFacade(GainCacheFacade&& other)  noexcept = default;
    GainCacheFacade & operator= (GainCacheFacade&& other)  noexcept = default;

    void assignQueryObjects(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) {
      _benefit_cache.assignQueryObjects(partIDs, pin_count_in_part, connectivity_sets);
      _penalty_cache.assignQueryObjects(partIDs, pin_count_in_part, connectivity_sets);
    }

    void parallel_resize(const HypernodeID num_nodes, const PartitionID k) {
        // make sure, this is called only if the gain cache was constructed with the empty constructor
        ASSERT(_num_nodes == 0 && _k == kInvalidPartition);

        _num_nodes = num_nodes;
        _k = k;
        tbb::parallel_invoke([&]{
            _benefit_cache.resize(std::string(GROUP_NAME),std::string(BENEFIT_KEY),num_nodes,k);
        }, [&] {
            _penalty_cache.resize(std::string(GROUP_NAME),std::string(PENALTY_KEY),num_nodes,k);
        });
    }

    size_t size_in_bytes() const {
        return _benefit_cache.size_in_bytes() + _penalty_cache.size_in_bytes();
    }

    size_t benefit_size() const {
        return _benefit_cache.size();
    }

    size_t penalty_size() const {
        return _penalty_cache.size();
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _benefit_cache.moveFromBenefit(u,m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, const PartitionID p, std::memory_order m = std::memory_order_relaxed) const {
        return _benefit_cache.moveFromBenefit(u,p,m);
    }

    // ! Not safe to call while moving/uncontracting
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefits, std::memory_order m = std::memory_order_relaxed) {
        _benefit_cache.storeRecomputedMoveFromBenefits(u, benefits, m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
        return _penalty_cache.moveToPenalty(u,to,m);
    }

    // ! Not safe to call while moving/uncontracting
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator,
                         const HyperedgeWeight incident_edges_weight, vec<HyperedgeWeight>& penalty_aggregator) {
        _penalty_cache.initializeEntry(u, incident_edges_weight, penalty_aggregator);
        _benefit_cache.initializeEntry(u, benefit_aggregator);
    }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> connectivity set queried on demand
    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v,
                                    const PartitionID block, const HypernodeID pin_count_in_part_after,
                                    IteratorRange<PinIteratorT> pins) {
        _penalty_cache.syncUpdateForUncontractCaseOne(he, we, v);
        _benefit_cache.syncUpdateForUncontractCaseOne(he, we, v, block, pin_count_in_part_after, pins);
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v, const PartitionID block,
                                    const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins,
                                    const ConnectivitySetSnapshot& connectivity_set, const PinCountInPart::Snapshot& pcip_snapshot) {
        _penalty_cache.asyncUpdateForUncontractCaseOne(we, v, connectivity_set);
        _benefit_cache.asyncUpdateForUncontractCaseOne(we, v, block, pin_count_in_part_after, pins, connectivity_set, pcip_snapshot);
    }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> connectivity set queried on demand
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u,
                                    const HypernodeID v) {
        _penalty_cache.syncUpdateForUncontractCaseTwo(he, we, u, v);
        _benefit_cache.syncUpdateForUncontractCaseTwo(he, we, u, v);
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v,
                                         const ConnectivitySetSnapshot& connectivity_set, const PinCountInPart::Snapshot& pcip_snapshot) {
        _penalty_cache.asyncUpdateForUncontractCaseTwo(we, u, v, connectivity_set);
        _benefit_cache.asyncUpdateForUncontractCaseTwo(we, u, v, connectivity_set, pcip_snapshot);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block) {
        _penalty_cache.updateForRestoringSinglePinNet(we, single_pin_in_he, block);
        _benefit_cache.updateForRestoringSinglePinNet(we, single_pin_in_he, block);
    }

    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
        _penalty_cache.updateForMove(we, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
        _benefit_cache.updateForMove(we, pins, from, pin_count_in_from_part_after, to, pin_count_in_to_part_after);
    }

private:

    static constexpr auto GROUP_NAME = "Refinement";
    static constexpr auto BENEFIT_KEY = "move_from_benefit";
    static constexpr auto PENALTY_KEY = "move_to_penalty";

    // ! Number of nodes in the cache
    HypernodeID _num_nodes;

    // ! Number of partitions in the cache
    PartitionID _k;

    // ! The cache for the move-from benefits
    MoveFromBenefitCache _benefit_cache;

    // ! The cache for the move-to penalties
    MoveToPenaltyCache _penalty_cache;

};

/// Move-from benefit cache where benefits for all blocks are aggregated, saving (k-1) * num_nodes memory compared to the
/// full benefit cache per node and block. This means that the benefit is generally invalid after one move and has to be
/// recomputed for the next move, though.
class AggregatedBenefitCache {

public:

    AggregatedBenefitCache() = default;

    AggregatedBenefitCache(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _part_ids(partIDs),
        _pins_in_parts(pin_count_in_part),
        _connectivity_sets(connectivity_sets),
        _move_from_benefit() {}

    AggregatedBenefitCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k,
                           const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) :
        _num_nodes(num_nodes),
        _k(k),
        _part_ids(partIDs),
        _pins_in_parts(pin_count_in_part),
        _connectivity_sets(connectivity_sets),
        _move_from_benefit(group, key, num_nodes, true, false) {}

    AggregatedBenefitCache(const AggregatedBenefitCache&) = delete;
    AggregatedBenefitCache & operator= (const AggregatedBenefitCache &) = delete;
    AggregatedBenefitCache(AggregatedBenefitCache&& other)  noexcept = default;
    AggregatedBenefitCache & operator= (AggregatedBenefitCache&& other)  noexcept = default;

    void assignQueryObjects(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) {
      _part_ids = partIDs;
      _pins_in_parts = pin_count_in_part;
      _connectivity_sets = connectivity_sets;
    }

    void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) {
        ASSERT(_num_nodes == 0 && _k == kInvalidPartition);
        _num_nodes = num_nodes;
        _k = k;
        _move_from_benefit.resize(group, key, num_nodes, true);
    }

    size_t size_in_bytes() const {
        return sizeof(CAtomic<HyperedgeWeight>) * _move_from_benefit.size();
    }

    size_t size() const {
        return _move_from_benefit.size();
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _move_from_benefit[u].load(m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID, const PartitionID, std::memory_order m = std::memory_order_relaxed) const {
        unused(m);
        ERROR("moveFromBenefit(HypernodeID, PartitionID) not supported on AggregatedBenefitCache.");
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator, std::memory_order m = std::memory_order_relaxed) {
      ASSERT(_part_ids);
        auto block_of_u = (*_part_ids)[u];
        _move_from_benefit[u].store(benefit_aggregator[block_of_u], m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator) {
      ASSERT(_part_ids);
        auto block_of_u = (*_part_ids)[u];
        _move_from_benefit[u].store(benefit_aggregator[block_of_u], std::memory_order_relaxed);
        // Reset all entries to zero
        benefit_aggregator.assign(_k,0);
    }

    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseOne(const HyperedgeID, const HyperedgeWeight we, const HypernodeID v,
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
    }
    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight, const HypernodeID, const PartitionID,
                                    const HypernodeID, IteratorRange<PinIteratorT>,
                                         const ConnectivitySetSnapshot&) {
        ERROR("asyncUpdateForUncontractCaseOne(HyperedgeWeight, HypernodeID, PartitionID, HypernodeID, "
              "IteratorRange<PinIteratorT>, const ConnectivitySetSnapshot&) not supported on AggregatedBenefitCache.");
    }


//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u, const HypernodeID v) {
        // Since u is no longer incident to hyperedge he its contribution for decreasing
        // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).
        ASSERT(_part_ids);
        ASSERT(_pins_in_parts);
        auto block = (*_part_ids)[u];
        if ( _pins_in_parts->pinCountInPart(he, block) == 1 ) {
            _move_from_benefit[u].sub_fetch(we, std::memory_order_relaxed);
            _move_from_benefit[v].add_fetch(we, std::memory_order_relaxed);
        }
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight, const HypernodeID, const HypernodeID, const ConnectivitySetSnapshot&) {
        ERROR("asyncUpdateForUncontractCaseTwo(HyperedgeWeight, HypernodeID, HypernodeID, HypernodeID, "
              "const ConnectivitySetSnapshot&) not supported on AggregatedBenefitCache.");
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID) {
        _move_from_benefit[single_pin_in_he].add_fetch(
                we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
      ASSERT(_part_ids);
        if (pin_count_in_from_part_after == 1) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, from);
                if ((*_part_ids)[u] == from) {
                    _move_from_benefit[u].fetch_add(we, std::memory_order_relaxed);
                }
            }
        }
        if (pin_count_in_to_part_after == 2) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, to);
                if ((*_part_ids)[u] == to) {
                    _move_from_benefit[u].fetch_sub(we, std::memory_order_relaxed);
                }
            }
        }
    }

private:

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
        unused(u);
        unused(p);
        ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
        ASSERT(p != kInvalidPartition && p < _k);
        ASSERT(isValidEntry(u,p));
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidEntry(const HypernodeID u, const PartitionID) const {
        return u < size();
    }

    // ! Total number of nodes that penalties are stored for
    HypernodeID _num_nodes;

    // ! Number of partitions that penalties are stored for
    PartitionID _k;

    const PartIDs* _part_ids;
    const PinCountInPart* _pins_in_parts;
    const ConnectivitySets* _connectivity_sets;

    // ! For each node and block, the sum of weights of incident edges with exactly one pin in that block
    Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

};

class FullBenefitCache {

public:

    FullBenefitCache() = default;

    FullBenefitCache(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets, parallel_tag_t) :
            _num_nodes(0),
            _k(kInvalidPartition),
            _part_ids(partIDs),
            _pins_in_parts(pin_count_in_part),
            _connectivity_sets(connectivity_sets),
            _move_from_benefit() {}

    FullBenefitCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k,
                     const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) :
            _num_nodes(num_nodes),
            _k(k),
            _part_ids(partIDs),
            _pins_in_parts(pin_count_in_part),
            _connectivity_sets(connectivity_sets),
            _move_from_benefit(group, key, size_t(num_nodes) * size_t(k), true, false) {}

    FullBenefitCache(const FullBenefitCache&) = delete;
    FullBenefitCache & operator= (const FullBenefitCache &) = delete;
    FullBenefitCache(FullBenefitCache&& other)  noexcept = default;
    FullBenefitCache & operator= (FullBenefitCache&& other)  noexcept = default;

    void assignQueryObjects(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) {
      _part_ids = partIDs;
      _pins_in_parts = pin_count_in_part;
      _connectivity_sets = connectivity_sets;
    }

    void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) {
        ASSERT(_num_nodes == 0 && _k == kInvalidPartition);
        _num_nodes = num_nodes;
        _k = k;
        _move_from_benefit.resize(group, key, size_t(num_nodes) * size_t(k), true);
    }

    size_t size_in_bytes() const {
        return sizeof(CAtomic<HyperedgeWeight>) * _move_from_benefit.size();
    }

    size_t size() const {
        return _move_from_benefit.size();
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
      ASSERT(_part_ids);
      return _move_from_benefit[benefit_index(u, (*_part_ids)[u])].load(m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, const PartitionID p, std::memory_order m = std::memory_order_relaxed) const {
        return _move_from_benefit[benefit_index(u, p)].load(m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator, std::memory_order m = std::memory_order_relaxed) {
        for (PartitionID i = 0; i < _k; ++i) {
            _move_from_benefit[benefit_index(u,i)].store(benefit_aggregator[i], m);
        }
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator) {
        for (PartitionID i = 0; i < _k; ++i) {
            _move_from_benefit[benefit_index(u,i)].store(benefit_aggregator[i], std::memory_order_relaxed);
            // Reset entry to zero so aggregator can be reused
            benefit_aggregator[i] = 0;
        }
    }

    // Variant for concurrent moves -> use snapshots of pins and parts_with_one_pin
    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v,
                                    const PartitionID block,
                                    const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins,
                                         const ConnectivitySetSnapshot& connectivity_set,
                                         const PinCountInPart::Snapshot& pcip_snapshot) {

        // Calculate and add contribution of he to the benefit of the uncontracted node v
        asyncAddContributionOfHEToNodeBenefit(v, we, connectivity_set, pcip_snapshot);

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

    // Variant without concurrent moves and uncontractions on pins of this HE -> use pin range directly from incidence
    // array and query connectivity set
    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v,
                                    const PartitionID block,
                                    const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins) {

      // Calculate and add contribution of he to the benefit of the uncontracted node v
      syncAddContributionOfHEToNodeBenefit(v, he, we);

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

    // Variant for concurrent moves -> use snapshots of parts_with_one_pin
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v, const ConnectivitySetSnapshot& connectivity_set, const PinCountInPart::Snapshot& pcip_snapshot) {
        // In this case, u is replaced by v in hyperedge he
        // => Pin counts of hyperedge he does not change
        // Since u is no longer incident to hyperedge he its contribution for decreasing
        // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).

        // Remove benefits contributed by this hyperedge from u as it no longer belongs to the hyperedge
        asyncRemoveContributionOfHEFromNodeBenefit(u, we, connectivity_set, pcip_snapshot);

        // Add benefits contributed by this hyperedge to the benefit of the uncontracted node v
        asyncAddContributionOfHEToNodeBenefit(v, we, connectivity_set, pcip_snapshot);
    }

    // Variant without concurrent moves and uncontractions on pins of this HE -> query parts with one pin on demand
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u, const HypernodeID v) {
      // Remove benefits contributed by this hyperedge from u as it no longer belongs to the hyperedge
      syncRemoveContributionOfHEFromNodeBenefit(u, he, we);

      // Add benefits contributed by this hyperedge to the benefit of the uncontracted node v
      syncAddContributionOfHEToNodeBenefit(v, he, we);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block) {
        _move_from_benefit[benefit_index(single_pin_in_he, block)].add_fetch(
                we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
        if (pin_count_in_from_part_after == 0) {
            for (const auto& u : pins) {
                nodeGainAssertions(u, from);
                _move_from_benefit[benefit_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
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
            }
        } else if (pin_count_in_to_part_after == 2) {
            for (const auto& u : pins) {
                nodeGainAssertions(u, to);
                _move_from_benefit[benefit_index(u, to)].fetch_sub(we, std::memory_order_relaxed);
            }
        }
    }

private:

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    size_t benefit_index(const HypernodeID u, const PartitionID p) const {
        return size_t(u) * size_t(_k)  + size_t(p);
    }

    // Variant without concurrent moves and uncontractions on this HE -> query parts with one pin on demand
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncAddContributionOfHEToNodeBenefit(const HypernodeID v, const HyperedgeID he, const HyperedgeWeight we) {
      ASSERT(_connectivity_sets);
      ASSERT(_pins_in_parts);
      for (const auto & p : _connectivity_sets->connectivitySet(he)) {
        if (_pins_in_parts->pinCountInPart(he, p) == 1) {
          _move_from_benefit[benefit_index(v, p)].add_fetch(we, std::memory_order_relaxed);
        }
      }
    }

    // Variant with concurrent moves and uncontractions on this HE -> require snapshot of parts with one pin in the edge
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncAddContributionOfHEToNodeBenefit(HypernodeID v, HyperedgeWeight we, const ConnectivitySetSnapshot& connectivity_set, const PinCountInPart::Snapshot& pcip_snapshot) {
        // Add benefit for all blocks that include exactly one pin of he
        for (const auto& p : connectivity_set) {
            if (pcip_snapshot.pinCountInPart(p) == 1) {
              _move_from_benefit[benefit_index(v, p)].add_fetch(we, std::memory_order_relaxed);
            }
        }
    }

    // Variant without concurrent moves and uncontractions on this HE -> query parts with one pin on demand
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncRemoveContributionOfHEFromNodeBenefit(const HypernodeID u, const HyperedgeID he, const HyperedgeWeight we) {
      ASSERT(_connectivity_sets);
      ASSERT(_pins_in_parts);
      for (const auto & p : _connectivity_sets->connectivitySet(he)) {
        if (_pins_in_parts->pinCountInPart(he, p) == 1) {
          _move_from_benefit[benefit_index(u, p)].sub_fetch(we, std::memory_order_relaxed);
        }
      }
    }

    // Variant with concurrent moves and uncontractions on this HE -> require snapshot of parts with one pin in the edge
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncRemoveContributionOfHEFromNodeBenefit(HypernodeID u, HyperedgeWeight we, const ConnectivitySetSnapshot& connectivity_set, const PinCountInPart::Snapshot& pcip_snapshot) {
        // Remove benefit for all blocks in the connectivity set of he
        for (const auto& p : connectivity_set) {
          if (pcip_snapshot.pinCountInPart(p) == 1) {
            _move_from_benefit[benefit_index(u, p)].sub_fetch(we, std::memory_order_relaxed);
          }
        }
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
        unused(u);
        unused(p);
        ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
//        ASSERT(_hg_query_funcs->is_node_enabled(u), "Hypernode" << u << "is disabled");
        ASSERT(p != kInvalidPartition && p < _k);
        ASSERT(isValidEntry(u,p));
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidEntry(const HypernodeID u, const PartitionID p) const {
        return benefit_index(u, p) < size();
    }

    // ! Total number of nodes that benefits are stored for
    HypernodeID _num_nodes;

    // ! Number of partitions that benefits are stored for
    PartitionID _k;

    const PartIDs* _part_ids;
    const PinCountInPart* _pins_in_parts;
    const ConnectivitySets* _connectivity_sets;

    // ! For each node and block, the sum of weights of incident edges with exactly one pin in that block
    Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

};

/// Move-to penalty cache where a penalty value is stored for moving a hypernode to a block for every combination of node
/// and block
class FullPenaltyCache {

public:

    FullPenaltyCache() = default;

    FullPenaltyCache(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _part_ids(partIDs),
        _pins_in_parts(pin_count_in_part),
        _connectivity_sets(connectivity_sets),
        _move_to_penalty() {}

    FullPenaltyCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k,
                     const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) :
        _num_nodes(num_nodes),
        _k(k),
        _part_ids(partIDs),
        _pins_in_parts(pin_count_in_part),
        _connectivity_sets(connectivity_sets),
        _move_to_penalty(group, key, size_t(num_nodes) * size_t(k + 1), true, false) {}

    FullPenaltyCache(const FullPenaltyCache&) = delete;
    FullPenaltyCache & operator= (const FullPenaltyCache &) = delete;
    FullPenaltyCache(FullPenaltyCache&& other)  noexcept = default;
    FullPenaltyCache & operator= (FullPenaltyCache&& other)  noexcept = default;


    void assignQueryObjects(const PartIDs* partIDs, const PinCountInPart* pin_count_in_part, const ConnectivitySets* connectivity_sets) {
      _part_ids = partIDs;
      _pins_in_parts = pin_count_in_part;
      _connectivity_sets = connectivity_sets;
    }

    void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) {
        ASSERT(_num_nodes == 0 && _k == kInvalidPartition);
        _num_nodes = num_nodes;
        _k = k;
        _move_to_penalty.resize(group, key, size_t(num_nodes) * size_t(k + 1), true);
    }

    HypernodeID size_in_bytes() const {
        return sizeof(CAtomic<HyperedgeWeight>) * _move_to_penalty.size();
    }

    size_t size() const {
        return _move_to_penalty.size();
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
        return _move_to_penalty[incident_net_weight_index(u)].load(m) -
               _move_to_penalty[penalty_index(u, to)].load(m);
    }

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, const HyperedgeWeight incident_edges_weight, vec<HyperedgeWeight>& penalty_aggregator) {
        _move_to_penalty[incident_net_weight_index(u)].store(
                incident_edges_weight, std::memory_order_relaxed);
        for (PartitionID p = 0; p < _k; ++p) {
            _move_to_penalty[penalty_index(u,p)].store(
                    penalty_aggregator[p], std::memory_order_relaxed);
            penalty_aggregator[p] = 0;
        }
    }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> query connectivity set on demand

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v) {
        ASSERT(_connectivity_sets);
        // For all blocks contained in the connectivity set of hyperedge he
        // we increase the move_to_penalty for vertex v by w(e) =>
        // move_to_penalty is then w(I(v)) - move_to_penalty(v, p) for a
        // block p
        _move_to_penalty[incident_net_weight_index(v)].add_fetch(
                we, std::memory_order_relaxed);
        for (const PartitionID block : _connectivity_sets->connectivitySet(he)) {
            _move_to_penalty[penalty_index(v, block)].add_fetch(
                    we, std::memory_order_relaxed);
        }
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> require snapshot of connectivity set
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v, const ConnectivitySetSnapshot& connectivity_set) {
        // For all blocks contained in the connectivity set of hyperedge he
        // we increase the move_to_penalty for vertex v by w(e) =>
        // move_to_penalty is then w(I(v)) - move_to_penalty(v, p) for a
        // block p
        _move_to_penalty[incident_net_weight_index(v)].add_fetch(
                we, std::memory_order_relaxed);
        for (const PartitionID block : connectivity_set) {
            _move_to_penalty[penalty_index(v, block)].add_fetch(
                    we, std::memory_order_relaxed);
        }
    }


    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> query connectivity set on demand

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void syncUpdateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u, const HypernodeID v) {
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

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> require snapshot of connectivity set
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void asyncUpdateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v, const ConnectivitySetSnapshot& connectivity_set) {
        // For all blocks contained in the connectivity set of hyperedge he
        // we decrease the the move_to_penalty for vertex u and increase it for
        // vertex v by w(e)
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

//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block_of_single_pin) {
        _move_to_penalty[incident_net_weight_index(single_pin_in_he)].add_fetch(
                we, std::memory_order_relaxed);
        _move_to_penalty[penalty_index(single_pin_in_he, block_of_single_pin)].add_fetch(
                we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
//    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
        if (pin_count_in_from_part_after == 0) {
            for (const auto& u : pins) {
                nodeGainAssertions(u, from);
                _move_to_penalty[penalty_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
            }
        }

        if (pin_count_in_to_part_after == 1) {
            for (const auto& u : pins) {
                nodeGainAssertions(u, to);
                _move_to_penalty[penalty_index(u, to)].fetch_add(we, std::memory_order_relaxed);
            }
        }
    }

private:

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
        unused(u);
        unused(p);
        ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
//        ASSERT(_hg_query_funcs->is_node_enabled(u), "Hypernode" << u << "is disabled");
        ASSERT(p != kInvalidPartition && p < _k);
        ASSERT(isValidEntry(u,p));
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidEntry(const HypernodeID u, const PartitionID p) const {
        return penalty_index(u, p) < size();
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    size_t incident_net_weight_index(const HypernodeID u) const {
        return size_t(u) * ( _k + 1 );
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    size_t penalty_index(const HypernodeID u, const PartitionID p) const {
        return size_t(u) * ( _k + 1 )  + p + 1;
    }

    // ! Total number of nodes that penalties are stored for
    HypernodeID _num_nodes;

    // ! Number of partitions that penalties are stored for
    PartitionID _k;

    const PartIDs* _part_ids;
    const PinCountInPart* _pins_in_parts;
    const ConnectivitySets* _connectivity_sets;

    // ! For each node and block, the sum of weights of incident edges with zero pins in that block
    Array< CAtomic<HyperedgeWeight> > _move_to_penalty;

};


using LightGainCache = GainCacheFacade<AggregatedBenefitCache, FullPenaltyCache>;
using HeavyGainCache = GainCacheFacade<FullBenefitCache, FullPenaltyCache>;

} // namespace mt_kahypar::ds
//
// Created by mlaupichler on 02.06.21.
//

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/async/iterable_bit_set.h"

namespace mt_kahypar::ds {

    using PartIDFunc = std::function<PartitionID(const HypernodeID)>;
    using ConnSetFunc = std::function<IteratorRange<ConnectivitySets::Iterator> (const HyperedgeID)>;
    using NodeEnabledFunc = std::function<bool(const HypernodeID)>;
    using PinCountInPartFunc = std::function<HypernodeID(const HyperedgeID, const PartitionID)>;
    struct HGQueryFunctions {
        // ! Function to query the current partition of a hypernode
        PartIDFunc part_id;

        // ! Function to query the current connectivity-set of an hyperedge
        ConnSetFunc connectivity_set;

        // ! Function to query whether a hypernode is currently enabled
        NodeEnabledFunc is_node_enabled;

        // ! Function to query the number of pins of a Hyperedge in a part
        PinCountInPartFunc pin_count_in_part;
    };

template<class MoveFromBenefitCache = Mandatory, class MoveToPenaltyCache = Mandatory>
class GainCacheFacade {

public:

    // Empty constructor to satisfy no-parameter constructor of PartitionedHypergraph
    GainCacheFacade() = default;

    // Default empty constructor to get size zero gain cache to allow subsequent parallel construction with parallel_resize()
    explicit GainCacheFacade(HGQueryFunctions* hg_query_funcs, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _benefit_cache(hg_query_funcs, parallel_tag_t()),
        _penalty_cache(hg_query_funcs, parallel_tag_t()) {}

    GainCacheFacade(const HypernodeID num_nodes, const PartitionID k, HGQueryFunctions* hg_query_funcs) :
        _num_nodes(num_nodes),
        _k(k),
        _benefit_cache(std::string(GROUP_NAME),std::string(BENEFIT_KEY),num_nodes,k, hg_query_funcs),
        _penalty_cache(std::string(GROUP_NAME),std::string(BENEFIT_KEY),num_nodes,k, hg_query_funcs) {}

    GainCacheFacade(const GainCacheFacade&) = delete;
    GainCacheFacade & operator= (const GainCacheFacade &) = delete;
    GainCacheFacade(GainCacheFacade&& other)  noexcept = default;
    GainCacheFacade & operator= (GainCacheFacade&& other)  noexcept = default;

    void assignHGQueryFunctions(HGQueryFunctions* hg_query_funcs) {
        _penalty_cache.assignHGQueryFunctions(hg_query_funcs);
        _benefit_cache.assignHGQueryFunctions(hg_query_funcs);
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

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _benefit_cache.moveFromBenefit(u,m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, const PartitionID p, std::memory_order m = std::memory_order_relaxed) const {
        return _benefit_cache.moveFromBenefit(u,p,m);
    }

    // ! Not safe to call while moving/uncontracting
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefits, std::memory_order m = std::memory_order_relaxed) {
        _benefit_cache.storeRecomputedMoveFromBenefits(u, benefits, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
        return _penalty_cache.moveToPenalty(u,to,m);
    }

    // ! Not safe to call while moving/uncontracting
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator,
                         const HyperedgeWeight incident_edges_weight, vec<HyperedgeWeight>& penalty_aggregator) {
        _penalty_cache.initializeEntry(u, incident_edges_weight, penalty_aggregator);
        _benefit_cache.initializeEntry(u, benefit_aggregator);
    }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> connectivity set queried on demand
    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v,
                                    const PartitionID block, const HypernodeID pin_count_in_part_after,
                                    IteratorRange<PinIteratorT> pins) {
        _penalty_cache.updateForUncontractCaseOne(he, we, v);
        _benefit_cache.updateForUncontractCaseOne(we, v, block, pin_count_in_part_after, pins);
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v, const PartitionID block,
                                    const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins,
                                    PartitionBitSet& connectivity_set, PartitionBitSet& parts_with_one_pin) {
        _penalty_cache.updateForUncontractCaseOne(we, v, connectivity_set);
        _benefit_cache.updateForUncontractCaseOne(we, v, block, pin_count_in_part_after, pins, parts_with_one_pin);
    }

    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> connectivity set queried on demand
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u,
                                    const HypernodeID v,
                                    const HypernodeID pin_count_in_part_after) {
        _penalty_cache.updateForUncontractCaseTwo(he, we, u, v);
        _benefit_cache.updateForUncontractCaseTwo(we, u, v, pin_count_in_part_after);
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> snapshot of connectivity set
    // ! and pin counts given
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v,
                                    PartitionBitSet& connectivity_set, PartitionBitSet& parts_with_one_pin) {
        _penalty_cache.updateForUncontractCaseTwo(we, u, v, connectivity_set);
        _benefit_cache.updateForUncontractCaseTwo(we, u, v, parts_with_one_pin);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block) {
        _penalty_cache.updateForRestoringSinglePinNet(we, single_pin_in_he, block);
        _benefit_cache.updateForRestoringSinglePinNet(we, single_pin_in_he, block);
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
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

    AggregatedBenefitCache(HGQueryFunctions* hg_query_funcs, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _hg_query_funcs(hg_query_funcs),
        _move_from_benefit() {}

    AggregatedBenefitCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k, HGQueryFunctions* hg_query_funcs) :
        _num_nodes(num_nodes),
        _k(k),
        _hg_query_funcs(hg_query_funcs),
        _move_from_benefit(group, key, num_nodes, true, false) {}

    AggregatedBenefitCache(const AggregatedBenefitCache&) = delete;
    AggregatedBenefitCache & operator= (const AggregatedBenefitCache &) = delete;
    AggregatedBenefitCache(AggregatedBenefitCache&& other)  noexcept = default;
    AggregatedBenefitCache & operator= (AggregatedBenefitCache&& other)  noexcept = default;

    void assignHGQueryFunctions(HGQueryFunctions* hg_query_funcs) {
        _hg_query_funcs = hg_query_funcs;
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

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _move_from_benefit[u].load(m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID, const PartitionID, std::memory_order m = std::memory_order_relaxed) const {
        unused(m);
        ERROR("moveFromBenefit(HypernodeID, PartitionID) not supported on AggregatedBenefitCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeRecomputedMoveFromBenefits(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator, std::memory_order m = std::memory_order_relaxed) {
        auto block_of_u = _hg_query_funcs->part_id(u);
        _move_from_benefit[u].store(benefit_aggregator[block_of_u], m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator) {
        auto block_of_u = _hg_query_funcs->part_id(u);
        _move_from_benefit[u].store(benefit_aggregator[block_of_u], std::memory_order_relaxed);
        // Reset all entries to zero
        benefit_aggregator.assign(_k,0);
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v,
                                    const PartitionID block,
                                    const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins) {

        if ( pin_count_in_part_after == 2 ) {
            // u might be replaced by an other vertex in the batch
            // => search for other pin of the corresponding block and
            // subtract edge weight.
            for ( auto it = pins.begin(); it != pins.end(); ++it ) {
                const HypernodeID& pin = *it;
                if ( pin != v && _hg_query_funcs->part_id(pin) == block ) {
                    _move_from_benefit[pin].sub_fetch(we, std::memory_order_relaxed);
                    break;
                }
            }
        }
    }
    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeWeight, const HypernodeID, const PartitionID,
                                    const HypernodeID, IteratorRange<PinIteratorT>,
                                    PartitionBitSet&) {
        ERROR("updateForUncontractCaseOne(HyperedgeWeight, HypernodeID, PartitionID, HypernodeID, "
              "IteratorRange<PinIteratorT>, PartitionBitSet&) not supported on AggregatedBenefitCache.");
    }


    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v,
                                    const HypernodeID pin_count_in_block) {
        // Since u is no longer incident to hyperedge he its contribution for decreasing
        // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).
        if ( pin_count_in_block == 1 ) {
            _move_from_benefit[u].sub_fetch(we, std::memory_order_relaxed);
            _move_from_benefit[v].add_fetch(we, std::memory_order_relaxed);
        }
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeWeight, const HypernodeID, const HypernodeID, PartitionBitSet&) {
        ERROR("updateForUncontractCaseTwo(HyperedgeWeight, HypernodeID, HypernodeID, HypernodeID, "
              "PartitionBitSet&) not supported on AggregatedBenefitCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID) {
        _move_from_benefit[single_pin_in_he].add_fetch(
                we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
        if (pin_count_in_from_part_after == 1) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, from);
                if (_hg_query_funcs->part_id(u) == from) {
                    _move_from_benefit[u].fetch_add(we, std::memory_order_relaxed);
                }
            }
        }
        if (pin_count_in_to_part_after == 2) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, to);
                if (_hg_query_funcs->part_id(u) == to) {
                    _move_from_benefit[u].fetch_sub(we, std::memory_order_relaxed);
                }
            }
        }
    }

private:

    template<typename PinIteratorT, typename F>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void doForEnabledPins(IteratorRange<PinIteratorT> pins, const F& f) {
        for (HypernodeID u : pins) {
            if (!_hg_query_funcs->is_node_enabled(u)) continue;
            f(u);
        }
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
        unused(u);
        unused(p);
        ASSERT(u < _num_nodes, "Hypernode" << u << "does not exist");
        ASSERT(_hg_query_funcs->is_node_enabled(u), "Hypernode" << u << "is disabled");
        ASSERT(p != kInvalidPartition && p < _k);
        ASSERT(isValidEntry(u,p));
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidEntry(const HypernodeID u, const PartitionID) const {
        return u < size();
    }

    // ! Total number of nodes that penalties are stored for
    HypernodeID _num_nodes;

    // ! Number of partitions that penalties are stored for
    PartitionID _k;

    // ! Functions used to query information about the current state of the hypergraph
    HGQueryFunctions* _hg_query_funcs;

    // ! For each node and block, the sum of weights of incident edges with exactly one pin in that block
    Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

};

class FullBenefitCache {

public:

    FullBenefitCache() = default;

    FullBenefitCache(HGQueryFunctions* hg_query_funcs, parallel_tag_t) :
            _num_nodes(0),
            _k(kInvalidPartition),
            _hg_query_funcs(hg_query_funcs),
            _move_from_benefit() {}

    FullBenefitCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k, HGQueryFunctions* hg_query_funcs) :
            _num_nodes(num_nodes),
            _k(k),
            _hg_query_funcs(hg_query_funcs),
            _move_from_benefit(group, key, size_t(num_nodes) * size_t(k), true, false) {}

    FullBenefitCache(const FullBenefitCache&) = delete;
    FullBenefitCache & operator= (const FullBenefitCache &) = delete;
    FullBenefitCache(FullBenefitCache&& other)  noexcept = default;
    FullBenefitCache & operator= (FullBenefitCache&& other)  noexcept = default;

    void assignHGQueryFunctions(HGQueryFunctions* hg_query_funcs) {
        _hg_query_funcs = hg_query_funcs;
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

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _move_from_benefit[benefit_index(u, _hg_query_funcs->part_id(u))].load(m);
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
    void initializeEntry(const HypernodeID u, vec<HyperedgeWeight>& benefit_aggregator) {
        for (PartitionID i = 0; i < _k; ++i) {
            _move_from_benefit[benefit_index(u,i)].store(benefit_aggregator[i], std::memory_order_relaxed);
            // Reset entry to zero so aggregator can be reused
            benefit_aggregator[i] = 0;
        }
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v,
                                    const PartitionID block,
                                    const HypernodeID pin_count_in_part_after, IteratorRange<PinIteratorT> pins,
                                    PartitionBitSet& parts_with_one_pin) {

        // Calculate and add contribution of he to the benefit of the uncontracted node v
        addContributionOfHEToNodeBenefit(v, we, parts_with_one_pin, true);

        // In this case, u and v are incident to hyperedge he after uncontraction
        if ( pin_count_in_part_after == 2 ) {
            // Reduce benefits for all pins in he for this block but exclude pins that are not enabled yet
            // (including v), as they never got the benefit to begin with
            for (const auto& pin : pins) {
                if (pin != v) {
                    _move_from_benefit[benefit_index(pin, block)].sub_fetch(we, std::memory_order_relaxed);
                }
            }
        }
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeWeight, const HypernodeID,
                                    const PartitionID,
                                    const HypernodeID, IteratorRange<PinIteratorT>) {

        ERROR("updateForUncontractCaseOne(HyperedgeWeight, HypernodeID, PartitionID, HypernodeID, "
              "IteratorRange<PinIteratorT> not supported on FullBenefitCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v, PartitionBitSet& parts_with_one_pin) {
        // In this case, u is replaced by v in hyperedge he
        // => Pin counts of hyperedge he does not change
        // Since u is no longer incident to hyperedge he its contribution for decreasing
        // the connectivity of he is shifted to vertex v => b(u) -= w(e), b(v) += w(e).

        // Remove benefits contributed by this hyperedge from u as it no longer belongs to the hyperedge
        removeContributionOfHEFromNodeBenefit(u, we, parts_with_one_pin, false);

        // Add benefits contributed by this hyperedge to the benefit of the uncontracted node v
        addContributionOfHEToNodeBenefit(v, we, parts_with_one_pin, true);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) {
        ERROR("updateForUncontractCaseTwo(HyperedgeWeight, HypernodeID, HypernodeID, HypernodeID, "
              "not supported on FullBenefitCache.");
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block) {
        _move_from_benefit[benefit_index(single_pin_in_he, block)].add_fetch(
                we, std::memory_order_relaxed);
    }

    template<class PinIteratorT>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForMove(const HyperedgeWeight we, IteratorRange<PinIteratorT> pins,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {
        if (pin_count_in_from_part_after == 0) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, from);
                _move_from_benefit[benefit_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
            }
        } else if (pin_count_in_from_part_after == 1) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, from);
                _move_from_benefit[benefit_index(u, from)].fetch_add(we, std::memory_order_relaxed);
            }
        }

        if (pin_count_in_to_part_after == 1) {
            for (const auto& u :pins) {
                nodeGainAssertions(u, to);
                _move_from_benefit[benefit_index(u, to)].fetch_add(we, std::memory_order_relaxed);
            }
        } else if (pin_count_in_to_part_after == 2) {
            for (const auto& u :pins) {
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

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addContributionOfHEToNodeBenefit(HypernodeID v, HyperedgeWeight we, PartitionBitSet& parts_with_one_pin,
                                          bool clear_bitset = false) {
        // Add benefit for all blocks that include exactly one pin of he
        for (const auto& p : parts_with_one_pin) {
            _move_from_benefit[benefit_index(v, p)].add_fetch(we, std::memory_order_relaxed);
            if (clear_bitset) parts_with_one_pin.set_false(p);
        }
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void removeContributionOfHEFromNodeBenefit(HypernodeID u, HyperedgeWeight we, PartitionBitSet& parts_with_one_pin,
                                               bool clear_bitset = false) {
        // Remove benefit for all blocks in the connectivity set of he
        for (const auto& p : parts_with_one_pin) {
            _move_from_benefit[benefit_index(u, p)].sub_fetch(we, std::memory_order_relaxed);
            if (clear_bitset) parts_with_one_pin.set_false(p);
        }
    }

    template<typename PinIteratorT, typename F>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void doForEnabledPins(IteratorRange<PinIteratorT> pins, const F& f) {
        for (HypernodeID u : pins) {
            if (!_hg_query_funcs->is_node_enabled(u)) continue;
            f(u);
        }
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

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidEntry(const HypernodeID u, const PartitionID p) const {
        return benefit_index(u, p) < size();
    }

    // ! Total number of nodes that benefits are stored for
    HypernodeID _num_nodes;

    // ! Number of partitions that benefits are stored for
    PartitionID _k;

    // ! Functions used to query information about the current state of the hypergraph
    HGQueryFunctions* _hg_query_funcs;

    // ! For each node and block, the sum of weights of incident edges with exactly one pin in that block
    Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

};

/// Move-to penalty cache where a penalty value is stored for moving a hypernode to a block for every combination of node
/// and block
class FullPenaltyCache {

public:

    FullPenaltyCache() = default;

    FullPenaltyCache(HGQueryFunctions* hg_query_funcs, parallel_tag_t) :
        _num_nodes(0),
        _k(kInvalidPartition),
        _hg_query_funcs(hg_query_funcs),
        _move_to_penalty() {}

    FullPenaltyCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k, HGQueryFunctions* hg_query_funcs) :
        _num_nodes(num_nodes),
        _k(k),
        _hg_query_funcs(hg_query_funcs),
        _move_to_penalty(group, key, size_t(num_nodes) * size_t(k + 1), true, false) {}

    FullPenaltyCache(const FullPenaltyCache&) = delete;
    FullPenaltyCache & operator= (const FullPenaltyCache &) = delete;
    FullPenaltyCache(FullPenaltyCache&& other)  noexcept = default;
    FullPenaltyCache & operator= (FullPenaltyCache&& other)  noexcept = default;

    void assignHGQueryFunctions(HGQueryFunctions* hg_query_funcs) {
        _hg_query_funcs = hg_query_funcs;
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

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
        return _move_to_penalty[incident_net_weight_index(u)].load(std::memory_order_relaxed) -
               _move_to_penalty[penalty_index(u, to)].load(std::memory_order_relaxed);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
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
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID v) {
        // For all blocks contained in the connectivity set of hyperedge he
        // we increase the the move_to_penalty for vertex v by w(e) =>
        // move_to_penalty is than w(I(v)) - move_to_penalty(v, p) for a
        // block p
        _move_to_penalty[incident_net_weight_index(v)].add_fetch(
                we, std::memory_order_relaxed);
        for (const PartitionID block : _hg_query_funcs->connectivity_set(he)) {
            _move_to_penalty[penalty_index(v, block)].add_fetch(
                    we, std::memory_order_relaxed);
        }
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> require snapshot of connectivity set
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseOne(const HyperedgeWeight we, const HypernodeID v, PartitionBitSet& connectivity_set) {
        // For all blocks contained in the connectivity set of hyperedge he
        // we increase the the move_to_penalty for vertex v by w(e) =>
        // move_to_penalty is than w(I(v)) - move_to_penalty(v, p) for a
        // block p
        _move_to_penalty[incident_net_weight_index(v)].add_fetch(
                we, std::memory_order_relaxed);
        for (const PartitionID block : connectivity_set) {
            _move_to_penalty[penalty_index(v, block)].add_fetch(
                    we, std::memory_order_relaxed);
            connectivity_set.set_false(block);
        }
    }


    // ! Variant for synchronous uncoarsening -> no moves during uncontraction -> query connectivity set on demand
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeID he, const HyperedgeWeight we, const HypernodeID u, const HypernodeID v) {
        // For all blocks contained in the connectivity set of hyperedge he
        // we decrease the the move_to_penalty for vertex u and increase it for
        // vertex v by w(e)
        _move_to_penalty[incident_net_weight_index(u)].sub_fetch(
                we, std::memory_order_relaxed);
        _move_to_penalty[incident_net_weight_index(v)].add_fetch(
                we, std::memory_order_relaxed);
        for (const PartitionID block : _hg_query_funcs->connectivity_set(he)) {
            _move_to_penalty[penalty_index(u, block)].sub_fetch(
                    we, std::memory_order_relaxed);
            _move_to_penalty[penalty_index(v, block)].add_fetch(
                    we, std::memory_order_relaxed);
        }
    }

    // ! Variant for asynchronous uncoarsening -> moves possible during uncontraction -> require snapshot of connectivity set
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForUncontractCaseTwo(const HyperedgeWeight we, const HypernodeID u, const HypernodeID v, PartitionBitSet& connectivity_set) {
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
            connectivity_set.set_false(block);
        }
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void updateForRestoringSinglePinNet(const HyperedgeWeight we, const HypernodeID single_pin_in_he, const PartitionID block_of_single_pin) {
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
            for (const auto& u :pins) {
                nodeGainAssertions(u, from);
                _move_to_penalty[penalty_index(u, from)].fetch_sub(we, std::memory_order_relaxed);
            }
        }

        if (pin_count_in_to_part_after == 1) {
            for (const auto& u :pins) {
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

    // ! Functions used to query information about the current state of the hypergraph
    HGQueryFunctions* _hg_query_funcs;

    // ! For each node and block, the sum of weights of incident edges with zero pins in that block
    Array< CAtomic<HyperedgeWeight> > _move_to_penalty;

};

using LightGainCache = GainCacheFacade<AggregatedBenefitCache,FullPenaltyCache>;
using HeavyGainCache = GainCacheFacade<FullBenefitCache, FullPenaltyCache>;

} // namespace mt_kahypar::ds
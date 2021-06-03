//
// Created by mlaupichler on 02.06.21.
//

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar::ds {

template<class MoveFromBenefitCache = Mandatory, class MoveToPenaltyCache = Mandatory>
class GainCacheFacade {

public:

    // Default empty constructor to get size zero gain cache to allow subsequent parallel construction with parallel_resize()
    explicit GainCacheFacade(const TaskGroupID task_group_id) : _num_nodes(0), _k(kInvalidPartition), _benefit_cache(task_group_id), _penalty_cache(task_group_id) {}

    GainCacheFacade(const HypernodeID num_nodes, const PartitionID k) :
        _num_nodes(num_nodes),
        _k(k),
        _benefit_cache(std::string(GROUP_NAME),std::string(BENEFIT_KEY),num_nodes,k),
        _penalty_cache(std::string(GROUP_NAME),std::string(BENEFIT_KEY),num_nodes,k) {}

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
    bool isValidEntry(const HypernodeID u, const PartitionID p) const {
        return _benefit_cache.isValidEntry(u,p) && _penalty_cache.isValidEntry(u,p);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _benefit_cache.moveFromBenefit(u,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addFetchMoveFromBenefit(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _benefit_cache.addFetchMoveFromBenefit(u,p,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void subFetchMoveFromBenefit(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _benefit_cache.subFetchMoveFromBenefit(u,p,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeMoveFromBenefit(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _benefit_cache.storeMoveFromBenefit(u,p,we,m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
        return _penalty_cache.moveToPenalty(u,to,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addFetchMoveToPenalty(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _penalty_cache.addFetchMoveToPenalty(u,p,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void subFetchMoveToPenalty(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _penalty_cache.subFetchMoveToPenalty(u,p,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeMoveToPenalty(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _penalty_cache.storeMoveToPenalty(u,p,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addFetchIncidentNetsWeight(const HypernodeID u, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _penalty_cache.addFetchIncidentNetsWeight(u,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void subFetchIncidentNetsWeight(const HypernodeID u, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _penalty_cache.subFetchIncidentNetsWeight(u,we,m);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeIncidentNetsWeight(const HypernodeID u, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _penalty_cache.storeIncidentNetsWeight(u,we,m);
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

    explicit AggregatedBenefitCache(const TaskGroupID) : _move_from_benefit() {}
    AggregatedBenefitCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID) :
        _move_from_benefit(group, key, num_nodes, true, false) {}

    void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID) {
        _move_from_benefit.resize(group, key, num_nodes, true);
    }

    size_t size_in_bytes() const {
        return sizeof(CAtomic<HyperedgeWeight>) * _move_from_benefit.size();
    }

    size_t size() const {
        return _move_from_benefit.size();
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    bool isValidEntry(const HypernodeID u, const PartitionID) const {
        return u < size();
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveFromBenefit(const HypernodeID u, std::memory_order m = std::memory_order_relaxed) const {
        return _move_from_benefit[u].load(m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addFetchMoveFromBenefit(const HypernodeID u, const PartitionID, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_from_benefit[u].add_fetch(we, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void subFetchMoveFromBenefit(const HypernodeID u, const PartitionID, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_from_benefit[u].sub_fetch(we, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeMoveFromBenefit(const HypernodeID u, const PartitionID, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_from_benefit[u].store(we,m);
    }

private:

    // ! For each node and block, the sum of weights of incident edges with exactly one pin in that block
    Array< CAtomic<HyperedgeWeight> > _move_from_benefit;

};

/// Move-to penalty cache where a penalty value is stored for moving a node to a block for every combination of node
/// and block
class FullPenaltyCache {

public:

    explicit FullPenaltyCache(const TaskGroupID) : _num_nodes(0), _k(kInvalidPartition), _move_to_penalty() {}
    FullPenaltyCache(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) :
        _num_nodes(num_nodes),
        _k(k),
        _move_to_penalty(group, key, size_t(num_nodes) * size_t(k + 1), true, false) {}

    void resize(const std::string& group, const std::string& key, const HypernodeID num_nodes, const PartitionID k) {
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
    bool isValidEntry(const HypernodeID u, const PartitionID p) const {
        return penalty_index(u, p) < size();
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID to, std::memory_order m = std::memory_order_relaxed) const {
        return _move_to_penalty[incident_net_weight_index(u)].load(std::memory_order_relaxed) -
               _move_to_penalty[penalty_index(u, to)].load(std::memory_order_relaxed);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addFetchMoveToPenalty(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_to_penalty[penalty_index(u,p)].add_fetch(we, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void subFetchMoveToPenalty(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_to_penalty[penalty_index(u,p)].sub_fetch(we, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeMoveToPenalty(const HypernodeID u, const PartitionID p, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_to_penalty[penalty_index(u,p)].store(we, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void addFetchIncidentNetsWeight(const HypernodeID u, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_to_penalty[incident_net_weight_index(u)].add_fetch(we, m);
    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void subFetchIncidentNetsWeight(const HypernodeID u, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_to_penalty[incident_net_weight_index(u)].sub_fetch(we, m);

    }

    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void storeIncidentNetsWeight(const HypernodeID u, const HyperedgeWeight we, std::memory_order m = std::memory_order_relaxed) {
        _move_to_penalty[incident_net_weight_index(u)].store(we, m);
    }

private:

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

    // ! For each node and block, the sum of weights of incident edges with zero pins in that block
    Array< CAtomic<HyperedgeWeight> > _move_to_penalty;

};

using LightGainCache = GainCacheFacade<AggregatedBenefitCache,FullPenaltyCache>;

} // namespace mt_kahypar::ds
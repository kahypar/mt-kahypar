/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <queue>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/km1/km1_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_gain_computation.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

namespace rebalancer {
struct RebalancingMove {
    HypernodeID hn;
    PartitionID to;
    float priority;
};
}; // namespace rebalancer
template <typename GraphAndGainTypes>
class DeterministicRebalancer final : public IRebalancer {
private:
    using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
    using GainCache = typename GraphAndGainTypes::GainCache;
    using GainComputation = typename GraphAndGainTypes::GainComputation;
    using RatingMap = typename GainComputation::RatingMap;
    using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

public:

    explicit DeterministicRebalancer(HypernodeID num_nodes, const Context& context) :
        _context(context),
        _max_part_weights(nullptr),
        _current_k(context.partition.k),
        _gain_computation(context),
        _num_imbalanced_parts(0),
        _num_valid_targets(0),
        _part_weights(_context.partition.k),
        _moves(num_nodes) { std::cout << "SO DETERMINSITC " << std::endl;}

    explicit DeterministicRebalancer(HypernodeID num_nodes, const Context& context, GainCache&) :
        DeterministicRebalancer(num_nodes, context) {}

    explicit DeterministicRebalancer(HypernodeID num_nodes, const Context& context, gain_cache_t gain_cache) :
        DeterministicRebalancer(num_nodes, context, GainCachePtr::cast<GainCache>(gain_cache)) {}

    DeterministicRebalancer(const DeterministicRebalancer&) = delete;
    DeterministicRebalancer(DeterministicRebalancer&&) = delete;

    DeterministicRebalancer& operator= (const DeterministicRebalancer&) = delete;
    DeterministicRebalancer& operator= (DeterministicRebalancer&&) = delete;

    bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
        const vec<HypernodeID>&,
        Metrics& best_metrics,
        double) final;

    void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final {}

    bool refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t&,
        const vec<HypernodeID>&,
        vec<vec<Move>>&,
        Metrics&,
        const double) override final {
        ERR("deterministic rebalancer can not be used for unconstrained refinement");
    }

    bool refineAndOutputMovesLinearImpl(mt_kahypar_partitioned_hypergraph_t&,
        const vec<HypernodeID>&,
        vec<Move>&,
        Metrics&,
        const double) override final {
        ERR("deterministic rebalancer can not be used for unconstrained refinement");
    }
private:

    void resizeDataStructuresForCurrentK() {
        // If the number of blocks changes, we resize data structures
        // (can happen during deep multilevel partitioning)
        if (_current_k != _context.partition.k) {
            _current_k = _context.partition.k;
            _gain_computation.changeNumberOfBlocks(_current_k);
            _part_weights = parallel::scalable_vector<AtomicWeight>(_context.partition.k);
        }
    }

    void initializeDataStructures(const PartitionedHypergraph& phg);

    void updateImbalance(const PartitionedHypergraph& hypergraph, bool read_weights_from_graph);

    bool mayMoveNode(PartitionID block, HypernodeWeight hn_weight) const {
        double allowed_weight = _part_weights[block].load(std::memory_order_relaxed)
            - _context.partition.perfect_balance_part_weights[block];
        allowed_weight *= _context.refinement.deterministic_refinement.jet.heavy_vertex_exclusion_factor;
        return hn_weight <= allowed_weight;
    }

    HypernodeWeight imbalance(PartitionID block) const {
        return _part_weights[block].load(std::memory_order_relaxed) - _max_part_weights[block];
    }

    HypernodeWeight deadzoneForPart(PartitionID block) const {
        const HypernodeWeight balanced = _context.partition.perfect_balance_part_weights[block];
        const HypernodeWeight max = _max_part_weights[block];
        return max - _context.refinement.deterministic_refinement.jet.relative_deadzone_size * (max - balanced);
    }

    bool isValidTarget(const PartitionedHypergraph& hypergraph,
        PartitionID block,
        HypernodeWeight hn_weight,
        bool use_precise_part_weights) const {
        const HypernodeWeight block_weight = use_precise_part_weights ?
            hypergraph.partWeight(block) : _part_weights[block].load(std::memory_order_relaxed);
        return (block_weight < deadzoneForPart(block)) &&
            block_weight + hn_weight <= _max_part_weights[block];
    }

    rebalancer::RebalancingMove computeGainAndTargetPart(const PartitionedHypergraph& hypergraph,
        const HypernodeID hn,
        bool non_adjacent_blocks,
        bool use_precise_part_weights);

    bool changeNodePart(PartitionedHypergraph& phg,
        const HypernodeID hn,
        const PartitionID from,
        const PartitionID to,
        bool ensure_balanced) {
        // it happens spuriously that from == to, not entirely sure why (possibly due to moving heavy nodes)
        if (from == to || to == kInvalidPartition) {
            return false;
        }

        // This function is passed as lambda to the changeNodePart function and used
        // to calculate the "real" delta of a move (in terms of the used objective function).
        auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
            _gain_computation.computeDeltaForHyperedge(sync_update);
        };

        HypernodeWeight max_weight = ensure_balanced ? _max_part_weights[to] : std::numeric_limits<HypernodeWeight>::max();
        bool success = false;
        success = phg.changeNodePart(hn, from, to, max_weight, [] {}, objective_delta);
        ASSERT(success || ensure_balanced);
        return success;
    }

    void weakRebalancingRound(PartitionedHypergraph& phg);

    const Context& _context;
    const HypernodeWeight* _max_part_weights;
    PartitionID _current_k;
    GainComputation _gain_computation;
    PartitionID _num_imbalanced_parts;
    PartitionID _num_valid_targets;
    parallel::scalable_vector<AtomicWeight> _part_weights;
    parallel::scalable_vector<rebalancer::RebalancingMove> _moves;
};

}  // namespace kahypar

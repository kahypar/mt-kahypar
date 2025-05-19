/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Robert Krause <robert.krause@student.kit.edu>
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"

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
    explicit DeterministicRebalancer(HypernodeID, const Context& context) :
        _context(context),
        _max_part_weights(nullptr),
        _current_k(context.partition.k),
        _gain_computation(context),
        _num_imbalanced_parts(0),
        _moves(context.partition.k),
        _move_weights(context.partition.k),
        _tmp_potential_moves(context.partition.k),
        _previous_block_imbalance(context.partition.k),
        _block_has_only_heavy_vertices(context.partition.k) {}

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
        throw UnsupportedOperationException("deterministic rebalancer can not be used for unconstrained refinement");
    }

    bool refineAndOutputMovesLinearImpl(mt_kahypar_partitioned_hypergraph_t&,
        const vec<HypernodeID>&,
        vec<Move>&,
        Metrics&,
        const double) override final {
        throw UnsupportedOperationException("deterministic rebalancer can not be used for unconstrained refinement");
    }

private:

    void resizeDataStructuresForCurrentK() {
        // If the number of blocks changes, we resize data structures
        // (can happen during deep multilevel partitioning)
        if (_current_k != _context.partition.k) {
            _current_k = _context.partition.k;
            _gain_computation.changeNumberOfBlocks(_current_k);
            _moves.resize(_current_k);
            _move_weights.resize(_current_k);
            _tmp_potential_moves.resize(_current_k);
            _previous_block_imbalance.resize(_current_k);
            _block_has_only_heavy_vertices.resize(_current_k);
        }
    }

    void initializeDataStructures(const PartitionedHypergraph& phg);

    void updateImbalance(const PartitionedHypergraph& hypergraph);

    // ! decides wether the node is allowed to be moved based on the heavy vertex excclusion hyperparameter
    bool mayMoveNode(const PartitionedHypergraph& phg, PartitionID part, HypernodeWeight hn_weight) const {
        double allowed_weight = phg.partWeight(part) - _context.partition.perfect_balance_part_weights[part];
        allowed_weight *= _context.refinement.rebalancing.det_heavy_vertex_exclusion_factor;
        return hn_weight <= allowed_weight || _block_has_only_heavy_vertices[part];
    }

    HypernodeWeight imbalance(const PartitionedHypergraph& phg, PartitionID part) const {
        return phg.partWeight(part) - _max_part_weights[part];
    }

    HypernodeWeight deadzoneForPart(PartitionID part) const {
        const HypernodeWeight balanced = _context.partition.perfect_balance_part_weights[part];
        const HypernodeWeight max = _max_part_weights[part];
        return max - _context.refinement.rebalancing.det_relative_deadzone_size * (max - balanced);
    }

    bool isValidTarget(const PartitionedHypergraph& hypergraph, PartitionID part, HypernodeWeight hn_weight) const {
        const HypernodeWeight part_weight = hypergraph.partWeight(part);
        return (part_weight < deadzoneForPart(part)) &&
            part_weight + hn_weight <= _max_part_weights[part];
    }

    rebalancer::RebalancingMove computeGainAndTargetPart(const PartitionedHypergraph& hypergraph,
                                                         const HypernodeID hn,
                                                         bool non_adjacent_blocks);

    bool changeNodePart(PartitionedHypergraph& phg,
                        const HypernodeID hn,
                        const PartitionID from,
                        const PartitionID to,
                        bool ensure_balanced) {

        // This function is passed as lambda to the changeNodePart function and used
        // to calculate the "real" delta of a move (in terms of the used objective function).
        auto objective_delta = [&](const SynchronizedEdgeUpdate& sync_update) {
            _gain_computation.computeDeltaForHyperedge(sync_update);
        };

        HypernodeWeight max_weight = ensure_balanced ? _max_part_weights[to] : std::numeric_limits<HypernodeWeight>::max();
        bool success = false;
        success = PartitionedHypergraph::is_graph ? phg.changeNodePartNoSync(hn, from, to, max_weight) : phg.changeNodePart(hn, from, to, max_weight, [] {}, objective_delta);
        ASSERT(success || ensure_balanced);
        return success;
    }

    void weakRebalancingRound(PartitionedHypergraph& phg);

    bool checkPreviouslyOverweightParts(const PartitionedHypergraph& phg) const {
        for (PartitionID i = 0; i < _current_k; ++i) {
            const auto partWeight = phg.partWeight(i);
            unused(partWeight);
            if (_previous_block_imbalance[i] > 0) {
                ASSERT(partWeight <= _max_part_weights[i] || _block_has_only_heavy_vertices[i], V(partWeight) << V(_max_part_weights[i]));
            }
        }
        return true;
    }

    const Context& _context;
    const HypernodeWeight* _max_part_weights;
    PartitionID _current_k;
    GainComputation _gain_computation;
    PartitionID _num_imbalanced_parts;
    parallel::scalable_vector<parallel::scalable_vector<rebalancer::RebalancingMove>> _moves;
    parallel::scalable_vector<parallel::scalable_vector<HypernodeWeight>> _move_weights;
    parallel::scalable_vector<ds::StreamingVector<rebalancer::RebalancingMove>> _tmp_potential_moves;
    parallel::scalable_vector<HypernodeWeight> _previous_block_imbalance;
    parallel::scalable_vector<uint8_t> _block_has_only_heavy_vertices;
};

}  // namespace kahypar

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

#include "external_tools/parlaylib/include/parlay/primitives.h"


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
    using AttributedGains = typename GraphAndGainTypes::AttributedGains;
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
        _num_valid_targets(0),
        _moves(context.partition.k),
        _move_weights(context.partition.k),
        tmp_potential_moves(context.partition.k) {}

    explicit DeterministicRebalancer(HypernodeID num_nodes, const Context& context, GainCache&) :
        DeterministicRebalancer(num_nodes, context) {}

    explicit DeterministicRebalancer(HypernodeID num_nodes, const Context& context, gain_cache_t gain_cache) :
        DeterministicRebalancer(num_nodes, context, GainCachePtr::cast<GainCache>(gain_cache)) {}

    DeterministicRebalancer(const DeterministicRebalancer&) = delete;
    DeterministicRebalancer(DeterministicRebalancer&&) = delete;

    DeterministicRebalancer& operator= (const DeterministicRebalancer&) = delete;
    DeterministicRebalancer& operator= (DeterministicRebalancer&&) = delete;
    bool jetRebalanceImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
        Metrics& best_metrics,
        bool run_until_balanced) {
        return refineInternal(hypergraph, best_metrics, run_until_balanced);
    }

    bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
        const vec<HypernodeID>&,
        Metrics& best_metrics,
        double) {
        return refineInternal(hypergraph, best_metrics, true);
    };

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
    bool refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph, Metrics& best_metrics, bool run_until_balanced);

    void resizeDataStructuresForCurrentK() {
        // If the number of blocks changes, we resize data structures
        // (can happen during deep multilevel partitioning)
        if (_current_k != _context.partition.k) {
            _current_k = _context.partition.k;
            _gain_computation.changeNumberOfBlocks(_current_k);
            _moves.resize(_current_k);
            _move_weights.resize(_current_k);
            tmp_potential_moves.resize(_current_k);
        }
    }

    void initializeDataStructures(const PartitionedHypergraph& phg);

    void updateImbalance(const PartitionedHypergraph& hypergraph);

    // ! decides wether the node is allowed to be moved based on the heavy vertex excclusion hyperparameter
    bool mayMoveNode(const PartitionedHypergraph& phg, PartitionID part, HypernodeWeight hn_weight) const {
        double allowed_weight = phg.partWeight(part) - _context.partition.perfect_balance_part_weights[part];
        allowed_weight *= _context.refinement.deterministic_refinement.jet.heavy_vertex_exclusion_factor;
        return hn_weight <= allowed_weight;
    }

    HypernodeWeight imbalance(const PartitionedHypergraph& phg, PartitionID part) const {
        return phg.partWeight(part) - _max_part_weights[part];
    }

    HypernodeWeight deadzoneForPart(PartitionID part) const {
        const HypernodeWeight balanced = _context.partition.perfect_balance_part_weights[part];
        const HypernodeWeight max = _max_part_weights[part];
        return max - _context.refinement.deterministic_refinement.jet.relative_deadzone_size * (max - balanced);
    }

    bool isValidTarget(const PartitionedHypergraph& hypergraph,
        PartitionID part,
        HypernodeWeight hn_weight) const {
        const HypernodeWeight part_weight = hypergraph.partWeight(part);
        return (part_weight < deadzoneForPart(part)) &&
            part_weight + hn_weight <= _max_part_weights[part];
    }

    rebalancer::RebalancingMove computeGainAndTargetPart(const PartitionedHypergraph& hypergraph,
        const HypernodeID hn,
        bool non_adjacent_blocks);

    template<bool isGraph>
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
        success = isGraph ? phg.changeNodePartNoSync(hn, from, to, max_weight) : phg.changeNodePart(hn, from, to, max_weight, [] {}, objective_delta);
        ASSERT(success || ensure_balanced);
        return success;
    }

    void weakRebalancingRound(PartitionedHypergraph& phg);

    void executeRebalancerForPartParallel(PartitionedHypergraph& phg, const size_t i, const size_t move_size, utils::Timer& timer) {
        // sort the moves from each overweight part by priority
        timer.start_timer("sorting", "Sorting");
        auto t0 = std::chrono::high_resolution_clock::now();
        parlay::sort_inplace(_moves[i], [&](const rebalancer::RebalancingMove& a, const rebalancer::RebalancingMove& b) {
            return a.priority < b.priority || (a.priority == b.priority && a.hn > b.hn);
        });
        auto t1 = std::chrono::high_resolution_clock::now();
        timer.stop_timer("sorting");
        t_sort += std::chrono::duration<double>
            (t1 - t0).count();
        // calculate perfix sum for each source-part to know which moves to execute (prefix_sum > current_weight - max_weight)
        timer.start_timer("find_moves", "Find Moves");
        t0 = std::chrono::high_resolution_clock::now();
        size_t last_move_idx = 0;
        if (move_size > _context.refinement.deterministic_refinement.jet.seq_find_rebalancing_moves) {
            if (move_size > _move_weights[i].size()) {
                _move_weights[i].resize(move_size);
            }
            tbb::parallel_for(0UL, move_size, [&](const size_t j) {
                _move_weights[i][j] = phg.nodeWeight(_moves[i][j].hn);
            });
            parallel_prefix_sum(_move_weights[i].begin(), _move_weights[i].begin() + move_size, _move_weights[i].begin(), std::plus<HypernodeWeight>(), 0);
            last_move_idx = std::upper_bound(_move_weights[i].begin(), _move_weights[i].begin() + move_size, phg.partWeight(i) - _max_part_weights[i] - 1) - _move_weights[i].begin();
        } else {
            HypernodeWeight sum = 0;
            for (; last_move_idx < move_size && sum <= phg.partWeight(i) - _max_part_weights[i] - 1; ++last_move_idx) {
                sum += phg.nodeWeight(_moves[i][last_move_idx].hn);
            }
            --last_move_idx;
        }
        t1 = std::chrono::high_resolution_clock::now();
        timer.stop_timer("find_moves");
        t_find += std::chrono::duration<double>
            (t1 - t0).count();

        timer.start_timer("exe_moves", "Execute Moves");
        t0 = std::chrono::high_resolution_clock::now();

        if (phg.is_graph) {
            tbb::parallel_for(0UL, last_move_idx + 1, [&](const size_t j) {
                const auto move = _moves[i][j];
                changeNodePart<true>(phg, move.hn, i, move.to, false);
            });
        } else {
            tbb::parallel_for(0UL, last_move_idx + 1, [&](const size_t j) {
                const auto move = _moves[i][j];
                changeNodePart<false>(phg, move.hn, i, move.to, false);
            });
        }
        t1 = std::chrono::high_resolution_clock::now();
        timer.stop_timer("exe_moves");
        t_exe += std::chrono::duration<double>
            (t1 - t0).count();
    }


    void executeRebalancerForPartSequential(PartitionedHypergraph& phg, const size_t i, const size_t move_size, utils::Timer& timer) {
        // sort the moves from each overweight part by priority
        timer.start_timer("sorting", "Sorting");
        auto t0 = std::chrono::high_resolution_clock::now();
        std::sort(_moves[i].begin(), _moves[i].begin() + move_size, [&](const rebalancer::RebalancingMove& a, const rebalancer::RebalancingMove& b) {
            return a.priority < b.priority || (a.priority == b.priority && a.hn > b.hn);
        });
        auto t1 = std::chrono::high_resolution_clock::now();
        timer.stop_timer("sorting");
        t_sort += std::chrono::duration<double>
            (t1 - t0).count();
        // calculate perfix sum for each source-part to know which moves to execute (prefix_sum > current_weight - max_weight)
        timer.start_timer("find_moves", "Find Moves");
        t0 = std::chrono::high_resolution_clock::now();
        HypernodeWeight sum = 0;
        size_t last_move_idx = 0;
        for (; last_move_idx < move_size && sum <= phg.partWeight(i) - _max_part_weights[i] - 1; ++last_move_idx) {
            sum += phg.nodeWeight(_moves[i][last_move_idx].hn);
        }
        t1 = std::chrono::high_resolution_clock::now();
        timer.stop_timer("find_moves");
        t_find += std::chrono::duration<double>
            (t1 - t0).count();

        timer.start_timer("exe_moves", "Execute Moves");
        t0 = std::chrono::high_resolution_clock::now();

        if (phg.is_graph) {
            for (size_t j = 0; j < last_move_idx; ++j) {
                const auto move = _moves[i][j];
                changeNodePart<true>(phg, move.hn, i, move.to, false);
            }
        } else {
            for (size_t j = 0; j < last_move_idx; ++j) {
                const auto move = _moves[i][j];
                changeNodePart<false>(phg, move.hn, i, move.to, false);
            }
        }
        t1 = std::chrono::high_resolution_clock::now();
        timer.stop_timer("exe_moves");
        t_exe += std::chrono::duration<double>
            (t1 - t0).count();
    }

    bool checkPreviouslyOverweightParts(const PartitionedHypergraph& phg)const {
        for (size_t i = 0; i < _moves.size(); ++i) {
            const auto partWeight = phg.partWeight(i);
            unused(partWeight);
            if (_moves[i].size() > 0) {
                ASSERT(partWeight >= deadzoneForPart(i) && partWeight <= _max_part_weights[i], V(partWeight) << V(deadzoneForPart(i)) << V(_max_part_weights[i]));
            }
        }
        return true;
    }

    const Context& _context;
    const HypernodeWeight* _max_part_weights;
    PartitionID _current_k;
    GainComputation _gain_computation;
    PartitionID _num_imbalanced_parts;
    PartitionID _num_valid_targets;
    parallel::scalable_vector<parallel::scalable_vector<rebalancer::RebalancingMove>> _moves;
    parallel::scalable_vector<parallel::scalable_vector<HypernodeWeight>> _move_weights;
    parallel::scalable_vector<ds::StreamingVector<rebalancer::RebalancingMove>> tmp_potential_moves;
    //parallel::scalable_vector<PartitionID> _part_before_round;

};

}  // namespace kahypar

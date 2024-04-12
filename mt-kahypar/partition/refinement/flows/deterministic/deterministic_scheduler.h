/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"
#include "mt-kahypar/partition/refinement/flows/problem_construction.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/participations_schedule.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_flow_refiner.h"
#include "mt-kahypar/partition/refinement/flows/flow_common.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_quotient_graph.h"

namespace mt_kahypar {

template<typename GraphAndGainTypes>
class DeterministicFlowRefinementScheduler final : public IRefiner {

    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

    using TypeTraits = typename GraphAndGainTypes::TypeTraits;
    using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
    using GainCache = typename GraphAndGainTypes::GainCache;
    using AttributedGains = typename GraphAndGainTypes::AttributedGains;

public:
    DeterministicFlowRefinementScheduler(const HypernodeID num_hypernodes,
        const HyperedgeID num_hyperedges,
        const Context& context,
        GainCache& gain_cache) :
        _context(context),
        _gain_cache(gain_cache),
        _scheduled_blocks(),
        _refiners(num_hypernodes, num_hyperedges, context),
        _quotient_graph(context),
        _schedule(context),
        _was_moved(num_hypernodes, uint8_t(false)) {}

    DeterministicFlowRefinementScheduler(const HypernodeID num_hypernodes,
        const HyperedgeID num_hyperedges,
        const Context& context,
        gain_cache_t gain_cache) :
        DeterministicFlowRefinementScheduler(num_hypernodes, num_hyperedges, context,
            GainCachePtr::cast<GainCache>(gain_cache)) {}

    DeterministicFlowRefinementScheduler(const DeterministicFlowRefinementScheduler&) = delete;
    DeterministicFlowRefinementScheduler(DeterministicFlowRefinementScheduler&&) = delete;

    DeterministicFlowRefinementScheduler& operator= (const DeterministicFlowRefinementScheduler&) = delete;
    DeterministicFlowRefinementScheduler& operator= (DeterministicFlowRefinementScheduler&&) = delete;


private:
    bool refineImpl(mt_kahypar_partitioned_hypergraph_t& phg,
        const vec<HypernodeID>& refinement_nodes,
        Metrics& metrics,
        double time_limit) final;

    void reportResults(const PartitionID i, const PartitionID j, const MoveSequence& sequence) {
        _schedule.reportResults(i, j, sequence);
    }

    HyperedgeWeight applyMoves(MoveSequence& sequence, PartitionedHypergraph& phg) {
        // Compute Part Weight Deltas
        vec<HypernodeWeight> part_weight_deltas(_context.partition.k, 0);
        for (Move& move : sequence.moves) {
            move.from = phg.partID(move.node);
            if (move.from != move.to) {
                const HypernodeWeight node_weight = phg.nodeWeight(move.node);
                part_weight_deltas[move.from] -= node_weight;
                part_weight_deltas[move.to] += node_weight;
            }
        }
        HyperedgeWeight improvement;
        vec<NewCutHyperedge> new_cut_hes;
        auto delta_func = [&](const SynchronizedEdgeUpdate& sync_update) {
            improvement -= AttributedGains::gain(sync_update);

            // Collect hyperedges with new blocks in its connectivity set
            if (sync_update.pin_count_in_to_part_after == 1) {
                // the corresponding block will be set in applyMoveSequence(...) function
                new_cut_hes.emplace_back(NewCutHyperedge{ sync_update.he, kInvalidPartition });
            }
        };
        applyMoveSequence(phg, _gain_cache, sequence, delta_func, _context.forceGainCacheUpdates(), _was_moved, new_cut_hes);
        assert(improvment == sequence.expected_improvement);
        addCutHyperedgesToQuotientGraph(new_cut_hes, phg);
        return improvement;
    }

    struct NewCutHyperedge {
        HyperedgeID he;
        PartitionID block;
    };

    template<typename F>
    void applyMoveSequence(PartitionedHypergraph& phg,
        GainCache& gain_cache,
        const MoveSequence& sequence,
        const F& objective_delta,
        const bool gain_cache_update,
        vec<uint8_t>& was_moved,
        vec<NewCutHyperedge>& new_cut_hes) {
        for (const Move& move : sequence.moves) {
            ASSERT(move.from == phg.partID(move.node));
            if (move.from != move.to) {
                changeNodePart(phg, gain_cache, move.node, move.from,
                    move.to, objective_delta, gain_cache_update);
                was_moved[move.node] = uint8_t(true);
                // If move increases the pin count of some hyperedges in block 'move.to' to one 1
                // we set the corresponding block here.
                int i = new_cut_hes.size() - 1;
                while (i >= 0 && new_cut_hes[i].block == kInvalidPartition) {
                    new_cut_hes[i].block = move.to;
                    --i;
                }
            }
        }
    }

    void addCutHyperedgesToQuotientGraph(const vec<NewCutHyperedge>& new_cut_hes, const PartitionedHypergraph& phg) {
        for (const NewCutHyperedge& new_cut_he : new_cut_hes) {
            ASSERT(new_cut_he.block != kInvalidPartition);
            _quotient_graph.addNewCutHyperedge(new_cut_he.he, new_cut_he.block, phg);
        }
    }

    template<typename F>
    bool changeNodePart(PartitionedHypergraph& phg,
        GainCache& gain_cache,
        const HypernodeID hn,
        const PartitionID from,
        const PartitionID to,
        const F& objective_delta,
        const bool gain_cache_update) {
        bool success = false;
        if (gain_cache_update && gain_cache.isInitialized()) {
            success = phg.changeNodePart(gain_cache, hn, from, to,
                std::numeric_limits<HypernodeWeight>::max(), [] {}, objective_delta);
        } else {
            success = phg.changeNodePart(hn, from, to,
                std::numeric_limits<HypernodeWeight>::max(), [] {}, objective_delta);
        }
        ASSERT(success);
        return success;
    }

    void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
        PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
        _quotient_graph.initialize(phg);
    }

    void resizeDataStructuresForCurrentK();

    const Context& _context;
    GainCache& _gain_cache;
    vec<BlockPair> _scheduled_blocks;
    tbb::enumerable_thread_specific<DeterministicFlowRefiner<GraphAndGainTypes>> _refiners;
    DeterministicQuotientGraph<TypeTraits> _quotient_graph;
    ParticipationsSchedule<TypeTraits> _schedule;
    vec<uint8_t> _was_moved;
};

}  // namespace kahypar

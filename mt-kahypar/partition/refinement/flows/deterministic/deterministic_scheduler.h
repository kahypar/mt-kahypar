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

#include "tbb/concurrent_queue.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"
#include "mt-kahypar/partition/refinement/flows/problem_construction.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/participations_schedule.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/participation_improvement_schedule.h"
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
        _refiners(),
        _quotient_graph(context),
        _schedule(context),
        _was_moved(num_hypernodes, uint8_t(false)),
        _apply_moves_lock(),
        num_hyperedges(num_hyperedges),
        num_hypernodes(num_hypernodes),
        _new_cut_hes(),
        _scheduled_blocks() {
        size_t numRefiners = std::min(size_t(_context.partition.k / 2), _context.shared_memory.num_threads);
        _refiners.reserve(numRefiners);
        for (size_t i = 0; i < numRefiners; ++i) {
            _refiners.emplace_back(std::make_unique<DeterministicFlowRefiner<GraphAndGainTypes>>(num_hypernodes, num_hyperedges, context));
        }
    }

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
        _apply_moves_lock.lock();
        // Compute Part Weight Deltas
        HyperedgeWeight improvement = 0;
        auto delta_func = [&](const SynchronizedEdgeUpdate& sync_update) {
            improvement -= AttributedGains::gain(sync_update);

            // Collect hyperedges with new blocks in its connectivity set
            if (sync_update.pin_count_in_to_part_after == 1) {
                // the corresponding block will be set in applyMoveSequence(...) function
                _new_cut_hes.emplace_back(NewCutHyperedge{ sync_update.he, kInvalidPartition });
            }
        };
        applyMoveSequence(phg, _gain_cache, sequence, delta_func, _context.forceGainCacheUpdates(), _was_moved);
        DBG << V(improvement) << ", " << V(sequence.expected_improvement);
        assert(improvement == sequence.expected_improvement);
        _apply_moves_lock.unlock();
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
        vec<uint8_t>& was_moved) {
        for (const Move& move : sequence.moves) {
            ASSERT(move.from == phg.partID(move.node));
            if (move.from != move.to) {
                changeNodePart(phg, gain_cache, move.node, move.from,
                    move.to, objective_delta, gain_cache_update);
                was_moved[move.node] = uint8_t(true);
                // If move increases the pin count of some hyperedges in block 'move.to' to one 1
                // we set the corresponding block here.
                int i = _new_cut_hes.size() - 1;
                while (i >= 0 && _new_cut_hes[i].block == kInvalidPartition) {
                    _new_cut_hes[i].block = move.to;
                    --i;
                }
            }
        }
    }

    void addCutHyperedgesToQuotientGraph(const PartitionedHypergraph& phg) {
        tbb::parallel_for(0UL, _new_cut_hes.size(), [&](const size_t i) {
            const NewCutHyperedge& new_cut_he = _new_cut_hes[i];
            _quotient_graph.addNewCutHyperedge(new_cut_he.he, new_cut_he.block, phg);
        });
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
    vec<std::unique_ptr<DeterministicFlowRefiner<GraphAndGainTypes>>> _refiners;
    DeterministicQuotientGraph<TypeTraits> _quotient_graph;
    ParticipationImprovementSchedule<TypeTraits> _schedule;
    vec<uint8_t> _was_moved;
    SpinLock _apply_moves_lock; // reset
    const HyperedgeID num_hyperedges;
    const HypernodeID num_hypernodes;
    vec<NewCutHyperedge> _new_cut_hes;
    tbb::concurrent_queue<ScheduledPair> _scheduled_blocks;
};

}  // namespace kahypar

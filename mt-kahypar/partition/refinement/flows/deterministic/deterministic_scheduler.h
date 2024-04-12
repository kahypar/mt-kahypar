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
        GainCache& gain_cache) : _scheduled_blocks(), _refiners(num_hypernodes, num_hyperedges, context), _quotient_graph(context), _schedule() {}

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

    void reportResults() {

    }

    void applyMoves() {

    }

    void initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
        PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
        _quotient_graph.initialize(phg);
    }

    void resizeDataStructuresForCurrentK();


    vec<BlockPair> _scheduled_blocks;
    tbb::enumerable_thread_specific<DeterministicFlowRefiner<GraphAndGainTypes>> _refiners;
    DeterministicQuotientGraph<TypeTraits> _quotient_graph;
    ParticipationsSchedule _schedule;
};

}  // namespace kahypar

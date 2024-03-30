/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/datastructures/buffered_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/reproducible_random.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_jet_refiner.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_label_propagation.h"

namespace mt_kahypar {

template<typename GraphAndGainTypes>
class CombinedDeterministicRefiner final : public IRefiner {

    using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
    using GainComputation = typename GraphAndGainTypes::GainComputation;
    using AttributedGains = typename GraphAndGainTypes::AttributedGains;
    using GainCache = typename GraphAndGainTypes::GainCache;
    using ActiveNodes = typename parallel::scalable_vector<HypernodeID>;
    using RatingMap = typename GainComputation::RatingMap;

public:

    explicit CombinedDeterministicRefiner(const HypernodeID num_hypernodes,
        const HyperedgeID num_hyperedges,
        const Context& context,
        gain_cache_t gain_cache,
        IRebalancer& rebalancer) : jetRefiner(num_hypernodes, num_hyperedges, context, gain_cache, rebalancer), lpRefiner(num_hypernodes, num_hyperedges, context, gain_cache, rebalancer) {}


private:
    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

    bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
        const vec<HypernodeID>& ignore,
        Metrics& best_metrics, double timeout) {
        lpRefiner.refine(hypergraph, ignore, best_metrics, timeout);
        jetRefiner.refine(hypergraph, ignore, best_metrics, timeout);
        return true;
    }

    void initializeImpl(mt_kahypar_partitioned_hypergraph_t& phg) {
        lpRefiner.initialize(phg);
        jetRefiner.initialize(phg);
    }

    DeterministicJetRefiner<GraphAndGainTypes> jetRefiner;
    DeterministicLabelPropagationRefiner<GraphAndGainTypes> lpRefiner;
};

}

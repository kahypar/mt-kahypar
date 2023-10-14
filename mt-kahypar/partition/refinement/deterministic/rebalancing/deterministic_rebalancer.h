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
template <typename GraphAndGainTypes>
class DeterministicRebalancer final : public IRebalancer {
private:
    using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
    using GainCache = typename GraphAndGainTypes::GainCache;
    using GainComputation = typename GraphAndGainTypes::GainComputation;
    using AtomicWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;

    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

public:

    explicit DeterministicRebalancer(const Context& context) :
        _context(context),
        _current_k(context.partition.k),
        _gain_computation(context),
        _part_weights(_context.partition.k) {}

    explicit DeterministicRebalancer(HypernodeID, const Context& context, GainCache&) :
        DeterministicRebalancer(context) {}

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

    const Context& _context;
    PartitionID _current_k;
    GainComputation _gain_computation;
    parallel::scalable_vector<AtomicWeight> _part_weights;
};

}  // namespace kahypar

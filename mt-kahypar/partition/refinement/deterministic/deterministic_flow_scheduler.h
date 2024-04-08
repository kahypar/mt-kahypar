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
#include "mt-kahypar/partition/refinement/flows/quotient_graph.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/utils/reproducible_random.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/flows/problem_construction.h"


namespace mt_kahypar {

template<typename GraphAndGainTypes>
class DeterministicFlowScheduler final : public IRefiner {

    using TypeTraits = typename GraphAndGainTypes::TypeTraits;
    using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
    using GainComputation = typename GraphAndGainTypes::GainComputation;
    using AttributedGains = typename GraphAndGainTypes::AttributedGains;
    using GainCache = typename GraphAndGainTypes::GainCache;

public:

    explicit DeterministicFlowScheduler(const HypernodeID num_hypernodes,
        const HyperedgeID num_hyperedges,
        const Context& context,
        GainCache& gain_cache) :
        _context(context),
        _gain_cache(gain_cache),
        _gain_computation(context, true /* disable_randomization */),
        _current_k(context.partition.k),
        _part_weights_lock(),
        _part_weights(context.partition.k, 0),
        _max_part_weights(context.partition.k, 0),
        _apply_moves_lock(),
        _quotient_graph(num_hyperedges, context),
        _constructor(num_hypernodes, num_hyperedges, context),
        _phg(nullptr),
        _was_moved(num_hypernodes, uint8_t(false)) {
        _refiner = FlowRefinementFactory::getInstance().createObject(FlowAlgorithm::flow_cutter, num_hyperedges, context);
    }

    DeterministicFlowScheduler(const HypernodeID num_hypernodes,
        const HyperedgeID num_hyperedges,
        const Context& context,
        gain_cache_t gain_cache) :
        DeterministicFlowScheduler(num_hypernodes, num_hyperedges, context,
            GainCachePtr::cast<GainCache>(gain_cache)) {}

    DeterministicFlowScheduler(const DeterministicFlowScheduler&) = delete;
    DeterministicFlowScheduler(DeterministicFlowScheduler&&) = delete;

    DeterministicFlowScheduler& operator= (const DeterministicFlowScheduler&) = delete;
    DeterministicFlowScheduler& operator= (DeterministicFlowScheduler&&) = delete;

    struct PartWeightUpdateResult {
        bool is_balanced = true;
        PartitionID overloaded_block = kInvalidPartition;
        HypernodeWeight overload_weight = 0;
    };

    /**
 * Applies the sequence of vertex moves to the partitioned hypergraph.
 * The method ensures that the move sequence does not violate
 * the balance constaint and not worsen solution quality.
 * Returns, improvement in solution quality.
 */
    HyperedgeWeight applyMoves(MoveSequence& sequence);

private:

    static constexpr bool debug = false;

    bool refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
        const vec<HypernodeID>& refinement_nodes,
        Metrics& best_metrics, double) final;

    void initializeImpl(mt_kahypar_partitioned_hypergraph_t&) final;

    void resizeDataStructuresForCurrentK();

    PartWeightUpdateResult partWeightUpdate(const vec<HypernodeWeight>& part_weight_deltas,
        const bool rollback);


    const Context& _context;
    GainCache& _gain_cache;
    GainComputation _gain_computation;
    PartitionID _current_k;

    // ! Maintains the part weights of each block
    SpinLock _part_weights_lock;
    vec<HypernodeWeight> _part_weights;
    vec<HypernodeWeight> _max_part_weights;

    SpinLock _apply_moves_lock;

    std::unique_ptr<IFlowRefiner> _refiner;

    // ! Contains information of all cut hyperedges between the
// ! blocks of the partition
    QuotientGraph<TypeTraits> _quotient_graph;

    // ! Responsible for construction of an flow problems
    ProblemConstruction<TypeTraits> _constructor;
    PartitionedHypergraph* _phg;

    // ! For each vertex it store wheather the corresponding vertex
// ! was moved or not
    vec<uint8_t> _was_moved;

};

}

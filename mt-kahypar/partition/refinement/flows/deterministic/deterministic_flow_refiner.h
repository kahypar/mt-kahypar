/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "algorithm/hyperflowcutter.h"
#include "algorithm/sequential_push_relabel.h"
#include "algorithm/parallel_push_relabel.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flows/sequential_construction.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_quotient_graph.h"
#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_problem_construction.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

template<typename GraphAndGainTypes>
class DeterministicFlowRefiner {

    static constexpr bool debug = false;
    static constexpr bool sequential = true;

    using PartitionedHypergraph = typename GraphAndGainTypes::PartitionedHypergraph;
    using TypeTraits = typename GraphAndGainTypes::TypeTraits;
public:
    explicit DeterministicFlowRefiner(const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,
        const Context& context) :
        _context(context),
        _flow_hg(),
        _sequential_hfc(_flow_hg, context.partition.seed),
        _parallel_hfc(_flow_hg, context.partition.seed),
        _sequential_construction(num_hyperedges, _flow_hg, _sequential_hfc, context),
        _problem_construction(num_hypernodes, num_hyperedges, context),
        _whfc_to_node() {
        _sequential_hfc.find_most_balanced = _context.refinement.flows.find_most_balanced_cut;
        _sequential_hfc.timer.active = false;
        _sequential_hfc.forceSequential(true);
        _sequential_hfc.setBulkPiercing(context.refinement.flows.pierce_in_bulk);

        _parallel_hfc.find_most_balanced = _context.refinement.flows.find_most_balanced_cut;
        _parallel_hfc.timer.active = false;
        _parallel_hfc.forceSequential(false);
        _parallel_hfc.setBulkPiercing(context.refinement.flows.pierce_in_bulk);
    }


    DeterministicFlowRefiner(const DeterministicFlowRefiner&) = delete;
    DeterministicFlowRefiner(DeterministicFlowRefiner&&) = delete;
    DeterministicFlowRefiner& operator= (const DeterministicFlowRefiner&) = delete;
    DeterministicFlowRefiner& operator= (DeterministicFlowRefiner&&) = delete;

    virtual ~DeterministicFlowRefiner() = default;

    void initialize(PartitionedHypergraph&) {
        _flow_hg.clear();
        _whfc_to_node.clear();
    }

    MoveSequence refine(PartitionedHypergraph& phg, DeterministicQuotientGraph<TypeTraits>& quotientGraph, const PartitionID block0, const PartitionID block1, const size_t seed) {
        return refineImpl(phg, quotientGraph, block0, block1, seed);
    }


private:

    MoveSequence refineImpl(PartitionedHypergraph& phg,
        DeterministicQuotientGraph<TypeTraits>& quotientGraph,
        const PartitionID block0,
        const PartitionID block1,
        const size_t seed) {
        MoveSequence sequence{ { }, 0 };
        Subhypergraph sub_hg = _problem_construction.construct(phg, quotientGraph, block0, block1);
        FlowProblem flow_problem = _sequential_construction.constructFlowHypergraph(phg, sub_hg, block0, block1, _whfc_to_node);
        bool flowcutter_succeeded = runFlowCutter(flow_problem, block0, block1, seed);
        if (flowcutter_succeeded) {
            extractMoveSequence(phg, flow_problem, sequence, block0, block1);
        }
        return sequence;

    }

    bool runFlowCutter(FlowProblem& flow_problem, const PartitionID block0, const PartitionID block1, const size_t seed) {
        whfc::Node s = flow_problem.source;
        whfc::Node t = flow_problem.sink;
        if (sequential) {
            _sequential_hfc.cs.setMaxBlockWeight(0, std::max(
                flow_problem.weight_of_block_0, _context.partition.max_part_weights[block0]));
            _sequential_hfc.cs.setMaxBlockWeight(1, std::max(
                flow_problem.weight_of_block_1, _context.partition.max_part_weights[block1]));

            _sequential_hfc.reset();
            _sequential_hfc.setSeed(seed);
            _sequential_hfc.setFlowBound(flow_problem.total_cut - flow_problem.non_removable_cut);
            return _sequential_hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t);
        } else {
            _parallel_hfc.cs.setMaxBlockWeight(0, std::max(
                flow_problem.weight_of_block_0, _context.partition.max_part_weights[block0]));
            _parallel_hfc.cs.setMaxBlockWeight(1, std::max(
                flow_problem.weight_of_block_1, _context.partition.max_part_weights[block1]));

            _parallel_hfc.reset();
            _parallel_hfc.setSeed(seed);
            _parallel_hfc.setFlowBound(flow_problem.total_cut - flow_problem.non_removable_cut);
            return _parallel_hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(s, t);
        }
    }

    void extractMoveSequence(const PartitionedHypergraph& phg, const FlowProblem& flow_problem, MoveSequence& sequence, const PartitionID block0, const PartitionID block1) {
        // We apply the solution if it either improves the cut or the balance of
        // the bipartition induced by the two blocks

        HyperedgeWeight new_cut = flow_problem.non_removable_cut;
        HypernodeWeight max_part_weight;
        if (sequential) {
            new_cut += _sequential_hfc.cs.flow_algo.flow_value;
            max_part_weight = std::max(_sequential_hfc.cs.source_weight, _sequential_hfc.cs.target_weight);
        } else {
            new_cut += _parallel_hfc.cs.flow_algo.flow_value;
            max_part_weight = std::max(_parallel_hfc.cs.source_weight, _parallel_hfc.cs.target_weight);
        }
        const bool improved_solution = new_cut < flow_problem.total_cut ||
            (new_cut == flow_problem.total_cut && max_part_weight < std::max(flow_problem.weight_of_block_0, flow_problem.weight_of_block_1));

        // Extract move sequence
        if (improved_solution) {
            sequence.expected_improvement = flow_problem.total_cut - new_cut;
            for (const whfc::Node& u : _flow_hg.nodeIDs()) {
                const HypernodeID hn = _whfc_to_node[u];
                if (hn != kInvalidHypernode) {
                    const PartitionID from = phg.partID(hn);
                    PartitionID to;
                    if (sequential) {
                        to = _sequential_hfc.cs.flow_algo.isSource(u) ? block0 : block1;
                    } else {
                        to = _parallel_hfc.cs.flow_algo.isSource(u) ? block0 : block1;
                    }

                    if (from != to) {
                        sequence.moves.push_back(Move{ from, to, hn, kInvalidGain });
                    }
                }
            }
        }
    }



    const Context& _context;
    FlowHypergraphBuilder _flow_hg;

    whfc::HyperFlowCutter<whfc::SequentialPushRelabel> _sequential_hfc;
    whfc::HyperFlowCutter<whfc::ParallelPushRelabel> _parallel_hfc;

    SequentialConstruction<GraphAndGainTypes> _sequential_construction;
    DeterministicProblemConstruction<TypeTraits> _problem_construction;
    //ParallelConstruction<GraphAndGainTypes> _parallel_construction;

    vec<HypernodeID> _whfc_to_node;

};
}  // namespace mt_kahypar

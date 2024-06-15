#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_scheduler.h"

namespace mt_kahypar {


template<typename GraphAndGainTypes>
bool DeterministicFlowRefinementScheduler<GraphAndGainTypes>::refineImpl(
    mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const parallel::scalable_vector<HypernodeID>&,
    Metrics& best_metrics,
    const double) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    if (num_hypernodes == phg.initialNumNodes()) {
        _schedule.setTopLevelFlag();
    }
    // if (!metrics::isBalanced(phg, _context)) {
    //     std::cout << "input is unbalanced" << std::endl;
    // } else {
    //     std::cout << "balanced" << std::endl;
    // }
    Metrics current_metrics = best_metrics;
    _schedule.initialize(hypergraph, _quotient_graph);
    HyperedgeWeight overall_delta = 0;
    HyperedgeWeight minImprovement = 0;
    std::atomic<HyperedgeWeight> round_delta(std::numeric_limits<HyperedgeWeight>::max());
    DBG << "";
    DBG << "------------------------------------------------------NEW LEVEL-----------------------------------------------------------------------";
    minImprovement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
    while (round_delta >= minImprovement && _schedule.hasActiveBlocks()) {
        //std::cout << V(round_delta) << ", " << V(minImprovement) << ", " << std::endl;
        size_t numScheduledBlocks = _schedule.getNextMatching(_scheduled_blocks, _quotient_graph);
        round_delta = 0;
        minImprovement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
        while (numScheduledBlocks > 0) {
            //vec<MoveSequence> sequences(_scheduled_blocks.size());
            tbb::parallel_for(0UL, _refiners.size(), [&](const size_t refinerIdx) {
                auto& refiner = *_refiners[refinerIdx];
                ScheduledPair sp;
                while (_scheduled_blocks.try_pop(sp)) {
                    //std::cout << "1" << std::endl;

                    // DeterministicFlowRefiner<GraphAndGainTypes> refiner(num_hypernodes, num_hyperedges,
                    //      _context);
                    refiner.initialize(phg);
                    //std::cout << "2" << std::endl;

                    MoveSequence moves = refiner.refine(phg, _quotient_graph, sp.bp.i, sp.bp.j, sp.seed);
                    //std::cout << "3" << std::endl;
                    //sequences[i] = moves;
                    const HyperedgeWeight improvement = applyMoves(moves, phg);
                    //std::cout << "4" << std::endl;

                    round_delta += improvement;
                    reportResults(sp.bp.i, sp.bp.j, moves);
                    _quotient_graph.reportImprovement(sp.bp.i, sp.bp.j, improvement);
                    //std::cout << "5" << std::endl;

                    //_solved_flow_problems++;
                    // for (auto v : moves.moves) {
                    //     std::cout << "(" << v.from << ", " << v.to << ", " << v.node << ", " << v.gain + ")";
                    // }
                    //std::cout << std::endl;
                }
            });
            //std::cout << "6" << std::endl;

            addCutHyperedgesToQuotientGraph(phg);
            //std::cout << "7" << std::endl;

            _new_cut_hes.clear();
            //tbb::parallel_for(0UL, sequences.size(), [&](const size_t i) {
            // for (size_t i = 0; i < sequences.size(); ++i) {
            //     //    for (size_t i = sequences.size() - 1; i < sequences.size(); --i) {
            //     const BlockPair& bp = _scheduled_blocks[i].bp;
            //     MoveSequence& moves = sequences[i];
            //     const HyperedgeWeight improvement = applyMoves(moves, phg);
            //     round_delta += improvement;
            //     reportResults(bp.i, bp.j, moves);
            //     _quotient_graph.reportImprovement(bp.i, bp.j, improvement);
            // }//);
            //_quotient_graph.initialize(phg);
            assert(metrics::isBalanced(phg, _context));
            //std::cout << V(_solved_flow_problems) << std::endl;
            DBG << "#################################################### NEXT MATCHING ######################################################";
            numScheduledBlocks = _schedule.getNextMatching(_scheduled_blocks, _quotient_graph);
        }
        overall_delta += round_delta;
        current_metrics.quality -= round_delta;
        DBG << "************************************************************ NEW ROUND *********************************************************";
        _schedule.resetForNewRound(_quotient_graph);
    }

    // Update metrics statistics
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality - overall_delta == metrics::quality(phg, _context),
        V(best_metrics.quality) << V(overall_delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality -= overall_delta;
    best_metrics.imbalance = metrics::imbalance(phg, _context);

    // Update Gain Cache
    if (_context.forceGainCacheUpdates() && _gain_cache.isInitialized()) {
        phg.doParallelForAllNodes([&](const HypernodeID& hn) {
            if (_was_moved[hn]) {
                _gain_cache.recomputeInvalidTerms(phg, hn);
                _was_moved[hn] = uint8_t(false);
            }
        });
    }
    return overall_delta > 0;
}

namespace {
#define DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER(X) DeterministicFlowRefinementScheduler<X>
}

INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER)

}
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
    utils::Timer& timer = utils::Utilities::instance().getTimer(_context.utility_id);

    Metrics current_metrics = best_metrics;
    timer.start_timer("scheduling_overhead", "Scheduling Overhead");
    _schedule.initialize(hypergraph, _quotient_graph);
    timer.stop_timer("scheduling_overhead");

    HyperedgeWeight overall_delta = 0;
    HyperedgeWeight minImprovement = 0;
    std::atomic<HyperedgeWeight> round_delta(std::numeric_limits<HyperedgeWeight>::max());
    DBG << "";
    DBG << "------------------------------------------------------NEW LEVEL-----------------------------------------------------------------------";
    minImprovement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
    while (round_delta >= minImprovement && _schedule.hasActiveBlocks()) {
        timer.start_timer("scheduling_overhead", "Scheduling Overhead");
        size_t numScheduledBlocks = _schedule.getNextMatching(_scheduled_blocks, _quotient_graph);
        timer.stop_timer("scheduling_overhead");
        round_delta = 0;
        minImprovement = _context.refinement.flows.min_relative_improvement_per_round * current_metrics.quality;
        while (numScheduledBlocks > 0) {
            timer.start_timer("flow_refiner", "Flow_Refiner");
            tbb::parallel_for(0UL, _refiners.size(), [&](const size_t refinerIdx) {
                auto& refiner = *_refiners[refinerIdx];
                ScheduledPair sp;
                while (_scheduled_blocks.try_pop(sp)) {
                    refiner.initialize(phg);

                    MoveSequence moves = refiner.refine(phg, _quotient_graph, sp.bp.i, sp.bp.j, sp.seed);
                    const HyperedgeWeight improvement = applyMoves(moves, phg);

                    round_delta += improvement;
                    reportResults(sp.bp.i, sp.bp.j, moves);
                    _quotient_graph.reportImprovement(sp.bp.i, sp.bp.j, improvement);

                }
            });
            timer.stop_timer("flow_refiner");
            addCutHyperedgesToQuotientGraph(phg);
            _new_cut_hes.clear();
            assert(metrics::isBalanced(phg, _context));
            DBG << "#################################################### NEXT MATCHING ######################################################";
            timer.start_timer("scheduling_overhead", "Scheduling Overhead");
            numScheduledBlocks = _schedule.getNextMatching(_scheduled_blocks, _quotient_graph);
            timer.stop_timer("scheduling_overhead");

        }
        overall_delta += round_delta;
        current_metrics.quality -= round_delta;
        DBG << "************************************************************ NEW ROUND *********************************************************";
        timer.start_timer("scheduling_overhead", "Scheduling Overhead");
        _schedule.resetForNewRound(_quotient_graph);
        timer.stop_timer("scheduling_overhead");
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
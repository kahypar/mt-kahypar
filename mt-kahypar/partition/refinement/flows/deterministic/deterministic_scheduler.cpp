#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_scheduler.h"

namespace mt_kahypar {


template<typename GraphAndGainTypes>
bool DeterministicFlowRefinementScheduler<GraphAndGainTypes>::refineImpl(
    mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const parallel::scalable_vector<HypernodeID>&,
    Metrics& best_metrics,
    const double) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    _schedule.initialize(hypergraph, _quotient_graph);
    size_t maxRounds = 5;
    std::atomic<HyperedgeWeight> overall_delta(0);
    for (size_t round = 0; round < maxRounds; ++round) {
        _scheduled_blocks = _schedule.getNextMatching(_quotient_graph);
        while (_scheduled_blocks.size() > 0) {
            tbb::parallel_for(0UL, _scheduled_blocks.size(), [&](const size_t i) {
                const BlockPair& bp = _scheduled_blocks[i];
                auto& refiner = _refiners.local();
                MoveSequence moves = refiner.refine(phg, _quotient_graph, bp.i, bp.j);
                const HyperedgeWeight improvement = applyMoves(moves, phg);
                overall_delta += improvement;
                reportResults(bp.i, bp.j, moves);
                _quotient_graph.reportImprovement(bp.i, bp.j, improvement);
            });
        }
        _schedule.resetForNewRound(_quotient_graph);
    }

    // Update metrics statistics
    HEAVY_REFINEMENT_ASSERT(best_metrics.quality + overall_delta == metrics::quality(phg, _context),
        V(best_metrics.quality) << V(overall_delta) << V(metrics::quality(phg, _context)));
    best_metrics.quality += overall_delta;
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
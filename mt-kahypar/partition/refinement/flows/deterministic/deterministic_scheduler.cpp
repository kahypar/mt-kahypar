#include "mt-kahypar/partition/refinement/flows/deterministic/deterministic_scheduler.h"

namespace mt_kahypar {


template<typename GraphAndGainTypes>
bool DeterministicFlowRefinementScheduler<GraphAndGainTypes>::refineImpl(
    mt_kahypar_partitioned_hypergraph_t& hypergraph,
    const parallel::scalable_vector<HypernodeID>&,
    Metrics& best_metrics,
    const double) {
    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    _schedule.initialize(hypergraph);
    _scheduled_blocks = _schedule.getNextRound();
    while (_scheduled_blocks.size() > 0) {
        tbb::parallel_for(0UL, _scheduled_blocks.size(), [&](const size_t i) {
            const BlockPair& bp = _scheduled_blocks[i];
            auto& refiner = _refiners.local();
            MoveSequence moves = refiner.refine(phg, _quotient_graph, bp.i, bp.j);
            reportResults();
            applyMoves();
        });
    }
}

namespace {
#define DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER(X) DeterministicFlowRefinementScheduler<X>
}

INSTANTIATE_CLASS_WITH_VALID_TRAITS(DETERMINISTIC_FLOW_REFINEMENT_SCHEDULER)

}
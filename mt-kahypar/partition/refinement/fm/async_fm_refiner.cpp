//
// Created by mlaupichler on 18.06.21.
//

#include <mt-kahypar/utils/stats.h>
#include "async_fm_refiner.h"

namespace mt_kahypar {

    template<typename FMStrategy>
    bool AsyncFMRefiner<FMStrategy>::refineImpl(PartitionedHypergraph &hypergraph,
                                                const parallel::scalable_vector <HypernodeID> &refinement_nodes,
                                                metrics::ThreadSafeMetrics &best_metrics, const double) {
      ASSERT(!refinement_nodes.empty(), "AsyncLPRefiner will not work without given seed refinement nodes. Cannot be used "
                                        "solely for rebalancing or for global refinement!");
      ASSERT(std::all_of(refinement_nodes.begin(),refinement_nodes.end(),[&](const HypernodeID& hn) {return hypergraph.nodeIsEnabled(hn);})
             && "Not all given seed nodes are enabled!");

      Gain delta = _fm->findMoves(hypergraph, refinement_nodes);

      if (debug && _context.type == kahypar::ContextType::main) {
        LOG << V(round) << V(delta) << V(metrics::km1(hypergraph)) << V(metrics::imbalance(hypergraph, _context))
            << V(refinement_nodes.size()) << _fm->stats.serialize();
      }

      // Update global part weight and sizes
      double imbalance;
      do {
        imbalance = metrics::imbalance(hypergraph, _context);
      } while (!best_metrics.update_imbalance_strong(imbalance));

      // Update metrics statistics
//      ASSERT(delta <= 0, "Async FM refiner worsen solution quality" << V(delta));
//      if (delta > 0) std::cout << "FM Delta > 0:\t" << delta << std::endl;
      best_metrics.fetch_add(delta, kahypar::Mode::direct_kway, kahypar::Objective::km1);
      utils::Stats::instance().update_stat("fm_improvement", std::abs(delta));
      return delta < 0;

    }

    template<typename FMStrategy>
    void AsyncFMRefiner<FMStrategy>::resetForGroup(ds::ContractionGroupID groupID) {
      _fm->resetForGroup(groupID);
    }

} // namespace mt_kahypar

#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h"

namespace mt_kahypar {
    template class AsyncFMRefiner<GainCacheStrategy<AsyncFMSharedData>>;
    template class AsyncFMRefiner<GainDeltaStrategy<AsyncFMSharedData>>;
    template class AsyncFMRefiner<RecomputeGainStrategy<AsyncFMSharedData>>;
    template class AsyncFMRefiner<GainCacheOnDemandStrategy<AsyncFMSharedData>>;
}
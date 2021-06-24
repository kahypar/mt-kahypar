//
// Created by mlaupichler on 18.06.21.
//

#pragma once

#include <mt-kahypar/partition/refinement/i_refiner.h>
#include <mt-kahypar/partition/refinement/async_refiners_common.h>
#include "async_kway_fm_core.h"

namespace mt_kahypar {


template<typename FMStrategy>
class AsyncFMRefiner : public IAsyncRefiner {

private:
    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = true;

public:

    AsyncFMRefiner(Hypergraph &hypergraph, const Context &context, ds::GroupLockManager *const lock_manager,
                   AsyncFMSharedData &shared_data) :
      _context(context),
      _fm(std::make_unique<AsyncKWayFM<FMStrategy>>(
          context, hypergraph.initialNumNodes(), shared_data, lock_manager)) {}

    int64_t getNumTotalAttemptedMoves() const override {
      ERROR("Number of moves currently not tracked in Async FM Refiner.");
    }

    int64_t getNumTotalMovedNodes() const override {
      ERROR("Number of moves currently not tracked in Async FM Refiner.");
    }

private:

    bool refineImpl(PartitionedHypergraph &hypergraph, const parallel::scalable_vector <HypernodeID> &refinement_nodes,
                    metrics::ThreadSafeMetrics &best_metrics, const double time_limit) override;

    void resetForGroup(ds::ContractionGroupID groupID) override;

    const Context& _context;
    std::unique_ptr<AsyncKWayFM<FMStrategy>> _fm;

};

} // namespace mt_kahypar




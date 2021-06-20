/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/


#include "kahypar/meta/registrar.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/do_nothing_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation/label_propagation_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation/async_lp_refiner.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_label_propagation.h"
#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"
#include "mt-kahypar/partition/refinement/fm/async_fm_refiner.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h"

#define REGISTER_LP_REFINER(id, refiner, t)                                                     \
  static kahypar::meta::Registrar<LabelPropagationFactory> JOIN(register_ ## refiner, t)(       \
    id,                                                                                         \
    [](Hypergraph& hypergraph, const Context& context) -> IRefiner* {                           \
    return new refiner(hypergraph, context);                                                    \
  })

#define REGISTER_ASYNC_LP_REFINER(id, refiner, t)                                                                       \
  static kahypar::meta::Registrar<AsyncLPRefinerFactory> JOIN(register_ ## refiner, t)(                                 \
    id,                                                                                                                 \
    [](Hypergraph& hypergraph, const Context& context,                                                                  \
    ds::GroupLockManager *lockManager, ds::ThreadSafeFlagArray<HypernodeID>* node_anti_duplicator,                      \
    ds::ThreadSafeFlagArray<HyperedgeID>* edge_anti_duplicator, AsyncNodeTracker *const fm_node_tracker)                \
    -> IAsyncRefiner* {                                                                                                 \
    return new refiner(hypergraph, context, lockManager, node_anti_duplicator, edge_anti_duplicator, fm_node_tracker);  \
  })

#define REGISTER_FM_REFINER(id, refiner, t)                                                     \
  static kahypar::meta::Registrar<FMFactory> JOIN(register_ ## refiner, t)(                     \
    id,                                                                                         \
    [](Hypergraph& hypergraph, const Context& context) -> IRefiner* {                           \
    return new refiner(hypergraph, context);                                                    \
  })

#define REGISTER_ASYNC_FM_REFINER(id, refiner, t)                                                 \
  static kahypar::meta::Registrar<AsyncFMRefinerFactory> JOIN(register_ ## refiner, t)(        \
      id,                                                                                         \
      [](Hypergraph &hypergraph, const Context &context, ds::GroupLockManager *const lock_manager,\
         AsyncFMSharedData &shared_data) -> IAsyncRefiner* {                                           \
        return new refiner(hypergraph, context, lock_manager, shared_data);                       \
      })

namespace mt_kahypar {
REGISTER_LP_REFINER(LabelPropagationAlgorithm::label_propagation_cut, LabelPropagationCutRefiner, Cut);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::label_propagation_km1, LabelPropagationKm1Refiner, Km1);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::deterministic, DeterministicLabelPropagationRefiner, Km1);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::do_nothing, DoNothingRefiner, 1);

using MultiTryKWayFMWithGainGache = MultiTryKWayFM<GainCacheStrategy<FMSharedData>>;
using MultiTryKWayFMWithGainGacheOnDemand = MultiTryKWayFM<GainCacheOnDemandStrategy<FMSharedData>>;
using MultiTryKWayFMWithGainDelta = MultiTryKWayFM<GainDeltaStrategy<FMSharedData>>;
using MultiTryKWayFMWithGainRecomputation = MultiTryKWayFM<RecomputeGainStrategy<FMSharedData>>;
REGISTER_FM_REFINER(FMAlgorithm::fm_gain_cache, MultiTryKWayFMWithGainGache, FMWithGainCache);
REGISTER_FM_REFINER(FMAlgorithm::fm_gain_cache_on_demand, MultiTryKWayFMWithGainGacheOnDemand, FMWithGainCacheOnDemand);
REGISTER_FM_REFINER(FMAlgorithm::fm_gain_delta, MultiTryKWayFMWithGainDelta, FMWithGainDelta);
REGISTER_FM_REFINER(FMAlgorithm::fm_recompute_gain, MultiTryKWayFMWithGainRecomputation, FMWithGainRecomputation);
REGISTER_FM_REFINER(FMAlgorithm::do_nothing, DoNothingRefiner, 2);

REGISTER_ASYNC_LP_REFINER(LabelPropagationAlgorithm::label_propagation_cut, AsyncLPCutRefiner, Cut);
REGISTER_ASYNC_LP_REFINER(LabelPropagationAlgorithm::label_propagation_km1, AsyncLPKm1Refiner, Km1);
REGISTER_ASYNC_LP_REFINER(LabelPropagationAlgorithm::do_nothing, DoNothingAsyncRefiner, 3);

using AsyncKWayFMWithGainCache = AsyncFMRefiner<GainCacheStrategy<AsyncFMSharedData>>;
using AsyncKWayFMWithGainGacheOnDemand = AsyncFMRefiner<GainCacheOnDemandStrategy<AsyncFMSharedData>>;
using AsyncKWayFMWithGainDelta = AsyncFMRefiner<GainDeltaStrategy<AsyncFMSharedData>>;
using AsyncKWayFMWithGainRecomputation = AsyncFMRefiner<RecomputeGainStrategy<AsyncFMSharedData>>;
REGISTER_ASYNC_FM_REFINER(FMAlgorithm::fm_gain_cache, AsyncKWayFMWithGainCache, FMWithGainCache);
REGISTER_ASYNC_FM_REFINER(FMAlgorithm::fm_gain_cache_on_demand, AsyncKWayFMWithGainGacheOnDemand, FMWithGainCacheOnDemand);
REGISTER_ASYNC_FM_REFINER(FMAlgorithm::fm_gain_delta, AsyncKWayFMWithGainDelta, FMWithGainDelta);
REGISTER_ASYNC_FM_REFINER(FMAlgorithm::fm_recompute_gain, AsyncKWayFMWithGainRecomputation, FMWithGainRecomputation);
REGISTER_ASYNC_FM_REFINER(FMAlgorithm::do_nothing, DoNothingAsyncRefiner, 4);

}  // namespace mt_kahypar

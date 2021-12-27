/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "kahypar/meta/registrar.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/do_nothing_refiner.h"
#include "mt-kahypar/partition/refinement/flows/do_nothing_refiner.h"
#include "mt-kahypar/partition/refinement/flows/flow_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation/label_propagation_refiner.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_label_propagation.h"
#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"
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

#define REGISTER_FM_REFINER(id, refiner, t)                                                     \
  static kahypar::meta::Registrar<FMFactory> JOIN(register_ ## refiner, t)(                     \
    id,                                                                                         \
    [](Hypergraph& hypergraph, const Context& context) -> IRefiner* {                           \
    return new refiner(hypergraph, context);                                                    \
  })

#define REGISTER_FLOW_REFINER(id, refiner, t)                                                 \
  static kahypar::meta::Registrar<FlowRefinementFactory> JOIN(register_ ## refiner, t)(       \
    id,                                                                                       \
    [](const Hypergraph& hypergraph, const Context& context) -> IFlowRefiner* {               \
    return new refiner(hypergraph, context);                                                  \
  })

namespace mt_kahypar {
REGISTER_LP_REFINER(LabelPropagationAlgorithm::label_propagation_cut, LabelPropagationCutRefiner, Cut);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::label_propagation_km1, LabelPropagationKm1Refiner, Km1);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::deterministic, DeterministicLabelPropagationRefiner, Km1);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::do_nothing, DoNothingRefiner, 1);

using MultiTryKWayFMWithGainGache = MultiTryKWayFM<GainCacheStrategy>;
using MultiTryKWayFMWithGainGacheOnDemand = MultiTryKWayFM<GainCacheOnDemandStrategy>;
using MultiTryKWayFMWithGainDelta = MultiTryKWayFM<GainDeltaStrategy>;
using MultiTryKWayFMWithGainRecomputation = MultiTryKWayFM<RecomputeGainStrategy>;
REGISTER_FM_REFINER(FMAlgorithm::fm_gain_cache, MultiTryKWayFMWithGainGache, FMWithGainCache);
REGISTER_FM_REFINER(FMAlgorithm::fm_gain_cache_on_demand, MultiTryKWayFMWithGainGacheOnDemand, FMWithGainCacheOnDemand);
REGISTER_FM_REFINER(FMAlgorithm::fm_gain_delta, MultiTryKWayFMWithGainDelta, FMWithGainDelta);
REGISTER_FM_REFINER(FMAlgorithm::fm_recompute_gain, MultiTryKWayFMWithGainRecomputation, FMWithGainRecomputation);
REGISTER_FM_REFINER(FMAlgorithm::do_nothing, DoNothingRefiner, 2);

REGISTER_FLOW_REFINER(FlowAlgorithm::do_nothing, DoNothingFlowRefiner, 3);
REGISTER_FLOW_REFINER(FlowAlgorithm::flow_cutter, FlowRefiner, Flows);

}  // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar/meta/abstract_factory.h"
#include "kahypar/meta/static_multi_dispatch_factory.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener.h"
#include "mt-kahypar/partition/coarsening/multilevel_coarsener.h"
#include "mt-kahypar/partition/coarsening/deterministic_multilevel_coarsener.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flows/i_flow_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation/label_propagation_refiner.h"
#include "mt-kahypar/partition/refinement/deterministic/deterministic_label_propagation.h"
#include "mt-kahypar/partition/refinement/fm/multitry_kway_fm.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h"
#include "mt-kahypar/partition/refinement/flows/scheduler.h"
#include "mt-kahypar/partition/refinement/flows/flow_refiner.h"

namespace mt_kahypar {

using CoarsenerFactory = kahypar::meta::Factory<CoarseningAlgorithm,
                                                ICoarsener* (*)(mt_kahypar_hypergraph_t, const Context&, uncoarsening_data_t*)>;

using MultilevelCoarsenerDispatcher = kahypar::meta::StaticMultiDispatchFactory<MultilevelCoarsener,
                                                                                ICoarsener,
                                                                                kahypar::meta::Typelist<TypeTraitsList,
                                                                                                        RatingScorePolicies,
                                                                                                        HeavyNodePenaltyPolicies,
                                                                                                        AcceptancePolicies> >;

using DeterministicCoarsenerDispatcher = kahypar::meta::StaticMultiDispatchFactory<DeterministicMultilevelCoarsener,
                                                                                   ICoarsener,
                                                                                   kahypar::meta::Typelist<TypeTraitsList>>;

using NLevelCoarsenerDispatcher = kahypar::meta::StaticMultiDispatchFactory<NLevelCoarsener,
                                                                            ICoarsener,
                                                                            kahypar::meta::Typelist<TypeTraitsList,
                                                                                                    RatingScorePolicies,
                                                                                                    HeavyNodePenaltyPolicies,
                                                                                                    AcceptancePolicies> >;

using LabelPropagationFactory = kahypar::meta::Factory<LabelPropagationAlgorithm,
                                                       IRefiner* (*)(HypernodeID, HyperedgeID, const Context&)>;

using Km1LabelPropagationDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                        LabelPropagationKm1Refiner,
                                        IRefiner,
                                        kahypar::meta::Typelist<TypeTraitsList>>;

using CutLabelPropagationDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                        LabelPropagationCutRefiner,
                                        IRefiner,
                                        kahypar::meta::Typelist<TypeTraitsList>>;

using DeterministicLabelPropagationDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                                  DeterministicLabelPropagationRefiner,
                                                  IRefiner,
                                                  kahypar::meta::Typelist<TypeTraitsList>>;

using FMFactory = kahypar::meta::Factory<FMAlgorithm,
                                         IRefiner* (*)(HypernodeID, HyperedgeID, const Context&)>;

template<typename TypeTraits>
using MultiTryKWayFMWithGainGache = MultiTryKWayFM<TypeTraits, GainCacheStrategy>;
using FMGainCacheDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                        MultiTryKWayFMWithGainGache,
                                        IRefiner,
                                        kahypar::meta::Typelist<TypeTraitsList>>;

template<typename TypeTraits>
using MultiTryKWayFMWithGainGacheOnDemand = MultiTryKWayFM<TypeTraits, GainCacheOnDemandStrategy>;
using FMGainCacheOnDemandDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                        MultiTryKWayFMWithGainGacheOnDemand,
                                        IRefiner,
                                        kahypar::meta::Typelist<TypeTraitsList>>;

template<typename TypeTraits>
using MultiTryKWayFMWithGainDelta = MultiTryKWayFM<TypeTraits, GainDeltaStrategy>;
using FMGainDeltaDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                        MultiTryKWayFMWithGainDelta,
                                        IRefiner,
                                        kahypar::meta::Typelist<TypeTraitsList>>;

template<typename TypeTraits>
using MultiTryKWayFMWithGainRecomputation = MultiTryKWayFM<TypeTraits, RecomputeGainStrategy>;
using FMGainRecomputationDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                        MultiTryKWayFMWithGainRecomputation,
                                        IRefiner,
                                        kahypar::meta::Typelist<TypeTraitsList>>;

using FlowSchedulerFactory = kahypar::meta::Factory<FlowAlgorithm,
                              IRefiner* (*)(const HypernodeID, const HyperedgeID, const Context&)>;
using FlowSchedulerDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                  FlowRefinementScheduler,
                                  IRefiner,
                                  kahypar::meta::Typelist<TypeTraitsList>>;

using FlowRefinementFactory = kahypar::meta::Factory<FlowAlgorithm,
                              IFlowRefiner* (*)(const HyperedgeID, const Context&)>;
using FlowRefinementDispatcher = kahypar::meta::StaticMultiDispatchFactory<
                                  FlowRefiner,
                                  IFlowRefiner,
                                  kahypar::meta::Typelist<TypeTraitsList>>;
}  // namespace mt_kahypar

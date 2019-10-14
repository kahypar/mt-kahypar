/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#pragma once

#include "kahypar/meta/abstract_factory.h"
#include "kahypar/meta/static_multi_dispatch_factory.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/preprocessing/i_community_assignment.h"
#include "mt-kahypar/partition/preprocessing/bin_packing_community_assignment.h"
#include "mt-kahypar/partition/preprocessing/policies/community_assignment_objective.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/community_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_community_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation_refiner.h"
#include "mt-kahypar/partition/refinement/policies/execution_policy.h"

namespace mt_kahypar {

using RedistributionFactory = kahypar::meta::Factory<CommunityAssignmentStrategy,
                                                     preprocessing::ICommunityAssignment* (*)(Hypergraph&, const Context&)>;

using BinPackingCommunityAssignmentDispatcher = kahypar::meta::StaticMultiDispatchFactory<preprocessing::BinPackingCommunityAssignment,
                                                                                          preprocessing::ICommunityAssignment,
                                                                                          kahypar::meta::Typelist<ObjectivePolicyClasses>>;

using CoarsenerFactory = kahypar::meta::Factory<CoarseningAlgorithm,
                                                ICoarsener* (*)(Hypergraph&, const Context&)>;

using CommunityCoarsenerDispatcher = kahypar::meta::StaticMultiDispatchFactory<CommunityCoarsener,
                                                                               ICoarsener,
                                                                               kahypar::meta::Typelist<RatingScorePolicies,
                                                                                                       HeavyNodePenaltyPolicies,
                                                                                                       AcceptancePolicies>>;


using LabelPropagationFactory = meta::Factory<LabelPropagationAlgorithm,
                                              IRefiner* (*)(Hypergraph&, const Context&)>;

using LabelPropagationKm1Dispatcher = kahypar::meta::StaticMultiDispatchFactory<LabelPropagationKm1Refiner,
                                                                                IRefiner,
                                                                                kahypar::meta::Typelist<ExecutionPolicyClasses>>;

using LabelPropagationCutDispatcher = kahypar::meta::StaticMultiDispatchFactory<LabelPropagationCutRefiner,
                                                                                IRefiner,
                                                                                kahypar::meta::Typelist<ExecutionPolicyClasses>>;




} // namespace mt_kahypar
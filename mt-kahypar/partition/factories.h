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

#include "tbb/task.h"

#include "kahypar/meta/abstract_factory.h"
#include "kahypar/meta/static_multi_dispatch_factory.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/partition/coarsening/community_coarsener.h"
#include "mt-kahypar/partition/coarsening/i_coarsener.h"
#include "mt-kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_community_policy.h"
#include "mt-kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/bin_packing_community_assignment.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/i_community_assignment.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/policies/community_assignment_objective.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation_refiner.h"
#include "mt-kahypar/partition/refinement/policies/execution_policy.h"

namespace mt_kahypar {
using RedistributionFactory = kahypar::meta::Factory<CommunityAssignmentStrategy,
                                                     preprocessing::ICommunityAssignment* (*)(Hypergraph&, const Context&)>;

using BinPackingCommunityAssignmentDispatcher = kahypar::meta::StaticMultiDispatchFactory<preprocessing::BinPackingCommunityAssignment,
                                                                                          preprocessing::ICommunityAssignment,
                                                                                          kahypar::meta::Typelist<ObjectivePolicyClasses> >;

using CoarsenerFactory = kahypar::meta::Factory<CoarseningAlgorithm,
                                                ICoarsener* (*)(Hypergraph&, const Context&, const TaskGroupID)>;

using CommunityCoarsenerDispatcher = kahypar::meta::StaticMultiDispatchFactory<CommunityCoarsener,
                                                                               ICoarsener,
                                                                               kahypar::meta::Typelist<RatingScorePolicies,
                                                                                                       HeavyNodePenaltyPolicies,
                                                                                                       AcceptancePolicies,
                                                                                                       ObjectivePolicyClasses> >;

using FlatInitialPartitionerFactory = kahypar::meta::Factory<InitialPartitioningAlgorithm,
                                                             tbb::task* (*)(tbb::task*, const InitialPartitioningAlgorithm, InitialPartitioningDataContainer&, const Context&)>;

using InitialPartitionerFactory = kahypar::meta::Factory<InitialPartitioningMode,
                                                         IInitialPartitioner* (*)(Hypergraph&, const Context&, const bool, const TaskGroupID)>;

using LabelPropagationFactory = kahypar::meta::Factory<LabelPropagationAlgorithm,
                                                       IRefiner* (*)(Hypergraph&, const Context&, const TaskGroupID)>;

using LabelPropagationKm1Dispatcher = kahypar::meta::StaticMultiDispatchFactory<LabelPropagationKm1Refiner,
                                                                                IRefiner,
                                                                                kahypar::meta::Typelist<ExecutionPolicyClasses> >;

using LabelPropagationCutDispatcher = kahypar::meta::StaticMultiDispatchFactory<LabelPropagationCutRefiner,
                                                                                IRefiner,
                                                                                kahypar::meta::Typelist<ExecutionPolicyClasses> >;
}  // namespace mt_kahypar

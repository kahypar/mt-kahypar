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

#pragma once

#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/preprocessing/sparsification/hypergraph_sparsifier.h"
#include "mt-kahypar/partition/preprocessing/sparsification/policies/similiar_net_combine.h"

#define REGISTER_DISPATCHED_COMMUNITY_ASSIGNER(id, dispatcher, ...)               \
  static kahypar::meta::Registrar<RedistributionFactory> register_ ## dispatcher( \
    id,                                                                           \
    [](Hypergraph& hypergraph, const Context& context) {                          \
    return dispatcher::create(                                                    \
      std::forward_as_tuple(hypergraph, context),                                 \
      __VA_ARGS__                                                                 \
      );                                                                          \
  })

#define REGISTER_HYPERGRAPH_SPARSIFIER(id, sparsifier)                                      \
  static kahypar::meta::Registrar<HypergraphSparsifierFactory> register_ ## sparsifier(     \
    id,                                                                                     \
    [](const Context& context, const TaskGroupID task_group_id)                             \
    -> IHypergraphSparsifier* {                                                             \
    return new sparsifier(context, task_group_id);                                          \
  })

namespace mt_kahypar {
REGISTER_DISPATCHED_COMMUNITY_ASSIGNER(CommunityAssignmentStrategy::bin_packing,
                                       BinPackingCommunityAssignmentDispatcher,
                                       kahypar::meta::PolicyRegistry<CommunityAssignmentObjective>::getInstance().getPolicy(
                                         context.preprocessing.community_redistribution.assignment_objective));

using HypergraphUnionSparsifier = HypergraphSparsifier<UnionCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::union_nets, HypergraphUnionSparsifier);
using HypergraphMaxSizeSparsifier = HypergraphSparsifier<MaxSizeCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::max_size, HypergraphMaxSizeSparsifier);
using HypergraphImportanceSparsifier = HypergraphSparsifier<NetImportanceCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::importance, HypergraphImportanceSparsifier);
using HypergraphUndefinedSparsifier = HypergraphSparsifier<UndefinedCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::UNDEFINED, HypergraphUndefinedSparsifier);
}  // namespace mt_kahypar

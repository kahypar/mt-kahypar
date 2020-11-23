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


#define REGISTER_DISPATCHED_COARSENER(id, dispatcher, ...)                                                      \
  static kahypar::meta::Registrar<CoarsenerFactory> register_ ## dispatcher(                                    \
    id,                                                                                                         \
    [](Hypergraph& hypergraph, const Context& context, const TaskGroupID task_group_id, const bool top_level) { \
    return dispatcher::create(                                                                                  \
      std::forward_as_tuple(hypergraph, context, task_group_id, top_level),                                     \
      __VA_ARGS__                                                                                               \
      );                                                                                                        \
  })

#define REGISTER_COARSENER(id, coarsener)                                                                                        \
  static kahypar::meta::Registrar<CoarsenerFactory> register_ ## coarsener(                                                      \
    id,                                                                                                                          \
    [](Hypergraph& hypergraph, const Context& context, const TaskGroupID task_group_id, const bool top_level) -> ICoarsener* {   \
    return new coarsener(hypergraph, context, task_group_id, top_level);                                                         \
  })


namespace mt_kahypar {
REGISTER_DISPATCHED_COARSENER(CoarseningAlgorithm::multilevel_coarsener,
                              MultilevelCoarsenerDispatcher,
                              kahypar::meta::PolicyRegistry<RatingFunction>::getInstance().getPolicy(
                                context.coarsening.rating.rating_function),
                              kahypar::meta::PolicyRegistry<HeavyNodePenaltyPolicy>::getInstance().getPolicy(
                                context.coarsening.rating.heavy_node_penalty_policy),
                              kahypar::meta::PolicyRegistry<AcceptancePolicy>::getInstance().getPolicy(
                                context.coarsening.rating.acceptance_policy));

REGISTER_DISPATCHED_COARSENER(CoarseningAlgorithm::nlevel_coarsener,
                              NLevelCoarsenerDispatcher,
                              kahypar::meta::PolicyRegistry<RatingFunction>::getInstance().getPolicy(
                                context.coarsening.rating.rating_function),
                              kahypar::meta::PolicyRegistry<HeavyNodePenaltyPolicy>::getInstance().getPolicy(
                                context.coarsening.rating.heavy_node_penalty_policy),
                              kahypar::meta::PolicyRegistry<AcceptancePolicy>::getInstance().getPolicy(
                                context.coarsening.rating.acceptance_policy));

}  // namespace mt_kahypar

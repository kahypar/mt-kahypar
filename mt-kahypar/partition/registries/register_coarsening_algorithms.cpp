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

#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/partition/coarsening/deterministic_multilevel_coarsener.h"


#define REGISTER_DISPATCHED_COARSENER(id, dispatcher, ...)                                      \
  static kahypar::meta::Registrar<CoarsenerFactory> register_ ## dispatcher(                    \
    id,                                                                                         \
    [](Hypergraph& hypergraph, const Context& context, UncoarseningData& uncoarseningData) {                  \
    return dispatcher::create(                                                                  \
      std::forward_as_tuple(hypergraph, context, uncoarseningData),                                    \
      __VA_ARGS__                                                                               \
      );                                                                                        \
  })

#define REGISTER_COARSENER(id, coarsener)                                                       \
  static kahypar::meta::Registrar<CoarsenerFactory> register_ ## coarsener(                     \
    id,                                                                                         \
    [](Hypergraph& hypergraph, const Context& context, UncoarseningData& uncoarseningData) -> ICoarsener* {   \
    return new coarsener(hypergraph, context, uncoarseningData);                                       \
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

REGISTER_COARSENER(CoarseningAlgorithm::deterministic_multilevel_coarsener, DeterministicMultilevelCoarsener);

}  // namespace mt_kahypar

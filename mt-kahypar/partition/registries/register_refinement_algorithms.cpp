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

#include "kahypar-resources/meta/registrar.h"
#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/do_nothing_refiner.h"
#include "mt-kahypar/partition/refinement/flows/do_nothing_refiner.h"

#define REGISTER_DISPATCHED_LP_REFINER(id, dispatcher, ...)                                            \
  static kahypar::meta::Registrar<LabelPropagationFactory> register_ ## dispatcher(                    \
    id,                                                                                                \
    [](const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,                             \
       const Context& context, gain_cache_t gain_cache) {                                              \
    return dispatcher::create(                                                                         \
      std::forward_as_tuple(num_hypernodes, num_hyperedges, context, gain_cache),                      \
      __VA_ARGS__                                                                                      \
      );                                                                                               \
  })

#define REGISTER_LP_REFINER(id, refiner, t)                                                      \
  static kahypar::meta::Registrar<LabelPropagationFactory> JOIN(register_ ## refiner, t)(        \
    id,                                                                                          \
    [](const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,                       \
       const Context& context, gain_cache_t gain_cache) -> IRefiner* {                           \
    return new refiner(num_hypernodes, num_hyperedges, context, gain_cache);                     \
  })

#define REGISTER_DISPATCHED_FM_REFINER(id, dispatcher, ...)                                            \
  static kahypar::meta::Registrar<FMFactory> register_ ## dispatcher(                                  \
    id,                                                                                                \
    [](const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,                             \
       const Context& context, gain_cache_t gain_cache) {                                              \
    return dispatcher::create(                                                                         \
      std::forward_as_tuple(num_hypernodes, num_hyperedges, context, gain_cache),                      \
      __VA_ARGS__                                                                                      \
      );                                                                                               \
  })

#define REGISTER_FM_REFINER(id, refiner, t)                                                      \
  static kahypar::meta::Registrar<FMFactory> JOIN(register_ ## refiner, t)(                      \
    id,                                                                                          \
    [](const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,                       \
       const Context& context, gain_cache_t gain_cache) -> IRefiner* {                           \
    return new refiner(num_hypernodes, num_hyperedges, context, gain_cache);                     \
  })

#define REGISTER_DISPATCHED_FLOW_SCHEDULER(id, dispatcher, ...)                                        \
  static kahypar::meta::Registrar<FlowSchedulerFactory> register_ ## dispatcher(                       \
    id,                                                                                                \
    [](const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,                             \
       const Context& context, gain_cache_t gain_cache) {                                              \
    return dispatcher::create(                                                                         \
      std::forward_as_tuple(num_hypernodes, num_hyperedges, context, gain_cache),                      \
      __VA_ARGS__                                                                                      \
      );                                                                                               \
  })

#define REGISTER_FLOW_SCHEDULER(id, refiner, t)                                                  \
  static kahypar::meta::Registrar<FlowSchedulerFactory> JOIN(register_ ## refiner, t)(           \
    id,                                                                                          \
    [](const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges,                       \
       const Context& context, gain_cache_t gain_cache) -> IRefiner* {                           \
    return new refiner(num_hypernodes, num_hyperedges, context, gain_cache);                     \
  })

#define REGISTER_DISPATCHED_REBALANCER(id, dispatcher, ...)                                            \
  static kahypar::meta::Registrar<RebalancerFactory> register_ ## dispatcher(                          \
    id,                                                                                                \
    [](const Context& context) {                                                                       \
    return dispatcher::create(                                                                         \
      std::forward_as_tuple(context),                                                                  \
      __VA_ARGS__                                                                                      \
      );                                                                                               \
  })

#define REGISTER_REBALANCER(id, refiner, t)                                                      \
  static kahypar::meta::Registrar<RebalancerFactory> JOIN(register_ ## refiner, t)(              \
    id,                                                                                          \
    [](const Context& context) -> IRefiner* {                                                    \
    return new refiner(context);                                                                 \
  })

#define REGISTER_DISPATCHED_FLOW_REFINER(id, dispatcher, ...)                                          \
  static kahypar::meta::Registrar<FlowRefinementFactory> register_ ## dispatcher(                      \
    id,                                                                                                \
    [](const HyperedgeID num_hyperedges, const Context& context) {                                     \
    return dispatcher::create(                                                                         \
      std::forward_as_tuple(num_hyperedges, context),                                                  \
      __VA_ARGS__                                                                                      \
      );                                                                                               \
  })

#define REGISTER_FLOW_REFINER(id, refiner, t)                                                   \
  static kahypar::meta::Registrar<FlowRefinementFactory> JOIN(register_ ## refiner, t)(         \
    id,                                                                                         \
    [](const HyperedgeID num_Hyperedges, const Context& context) -> IFlowRefiner* {             \
    return new refiner(num_Hyperedges, context);                                                \
  })

namespace mt_kahypar {
REGISTER_DISPATCHED_LP_REFINER(LabelPropagationAlgorithm::label_propagation,
                               LabelPropagationDispatcher,
                               kahypar::meta::PolicyRegistry<mt_kahypar_partition_type_t>::getInstance().getPolicy(
                                context.partition.partition_type),
                               kahypar::meta::PolicyRegistry<GainPolicy>::getInstance().getPolicy(
                                context.partition.gain_policy));
REGISTER_DISPATCHED_LP_REFINER(LabelPropagationAlgorithm::deterministic,
                               DeterministicLabelPropagationDispatcher,
                               kahypar::meta::PolicyRegistry<mt_kahypar_partition_type_t>::getInstance().getPolicy(
                                context.partition.partition_type));
REGISTER_LP_REFINER(LabelPropagationAlgorithm::do_nothing, DoNothingRefiner, 1);

REGISTER_DISPATCHED_FM_REFINER(FMAlgorithm::kway_fm,
                               FMDispatcher,
                               kahypar::meta::PolicyRegistry<mt_kahypar_partition_type_t>::getInstance().getPolicy(
                                context.partition.partition_type),
                               kahypar::meta::PolicyRegistry<GainPolicy>::getInstance().getPolicy(
                                context.partition.gain_policy));
REGISTER_FM_REFINER(FMAlgorithm::do_nothing, DoNothingRefiner, 2);

REGISTER_DISPATCHED_FLOW_SCHEDULER(FlowAlgorithm::flow_cutter,
                                   FlowSchedulerDispatcher,
                                   kahypar::meta::PolicyRegistry<mt_kahypar_partition_type_t>::getInstance().getPolicy(
                                    context.partition.partition_type),
                                   kahypar::meta::PolicyRegistry<GainPolicy>::getInstance().getPolicy(
                                     context.partition.gain_policy));
REGISTER_FLOW_SCHEDULER(FlowAlgorithm::do_nothing, DoNothingRefiner, 3);

REGISTER_DISPATCHED_REBALANCER(RebalancingAlgorithm::simple_rebalancer,
                               RebalancerDispatcher,
                               kahypar::meta::PolicyRegistry<mt_kahypar_partition_type_t>::getInstance().getPolicy(
                                context.partition.partition_type),
                               kahypar::meta::PolicyRegistry<GainPolicy>::getInstance().getPolicy(
                                context.partition.gain_policy));
REGISTER_REBALANCER(RebalancingAlgorithm::do_nothing, DoNothingRefiner, 4);

REGISTER_DISPATCHED_FLOW_REFINER(FlowAlgorithm::flow_cutter,
                                  FlowRefinementDispatcher,
                                  kahypar::meta::PolicyRegistry<mt_kahypar_partition_type_t>::getInstance().getPolicy(
                                   context.partition.partition_type),
                                  kahypar::meta::PolicyRegistry<GainPolicy>::getInstance().getPolicy(
                                    context.partition.gain_policy));
REGISTER_FLOW_REFINER(FlowAlgorithm::do_nothing, DoNothingFlowRefiner, 5);
}  // namespace mt_kahypar

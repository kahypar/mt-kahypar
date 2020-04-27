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
#include "mt-kahypar/partition/refinement/do_nothing_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation/label_propagation_refiner.h"

#define REGISTER_LP_REFINER(id, refiner, t)                                                                          \
  static kahypar::meta::Registrar<LabelPropagationFactory> JOIN(register_ ## refiner, t)(                            \
    id,                                                                                                              \
    [](PartitionedHypergraph<>& hypergraph, const Context& context, const TaskGroupID task_group_id) -> IRefiner* {  \
    return new refiner(hypergraph, context, task_group_id);                                                          \
  })

namespace mt_kahypar {
REGISTER_LP_REFINER(LabelPropagationAlgorithm::label_propagation_cut, LabelPropagationCutRefiner, Cut);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::label_propagation_km1, LabelPropagationKm1Refiner, Km1);
REGISTER_LP_REFINER(LabelPropagationAlgorithm::do_nothing, DoNothingRefiner, 1);
}  // namespace mt_kahypar

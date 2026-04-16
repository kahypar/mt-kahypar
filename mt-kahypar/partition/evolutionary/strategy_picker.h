/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2017 Robin Andre <robinandre1995@web.de>
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
#include <kahypar-resources/utils/randomize.h>

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/evo_partitioner.h"


namespace mt_kahypar::pick {
EvoMutateStrategy decideNextMutation(const Context& context, std::mt19937* rng = nullptr);
EvoDecision decideNextMove(const Context& context, std::mt19937* rng = nullptr) ;
ContextModifierParameters decideContextModificationParameters(const Context& context, std::mt19937* rng = nullptr);
} // namespace mt_kahypar::pick


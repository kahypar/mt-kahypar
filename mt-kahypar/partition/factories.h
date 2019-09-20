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
#include "mt-kahypar/partition/preprocessing/i_redistribution.h"
#include "mt-kahypar/partition/preprocessing/bin_packing_restribution.h"
#include "mt-kahypar/partition/preprocessing/policies/community_assignment_objective.h"

namespace mt_kahypar {

using RedistributionFactory = kahypar::meta::Factory<CommunityAssignmentStrategy,
                                                     preprocessing::IRedistribution* (*)(Hypergraph&, const Context&)>;

using BinPackingRedistributionDispatcher = kahypar::meta::StaticMultiDispatchFactory<preprocessing::BinPackingRedistribution,
                                                                                     preprocessing::IRedistribution,
                                                                                     kahypar::meta::Typelist<ObjectivePolicyClasses>>;

} // namespace mt_kahypar
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

#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

namespace mt_kahypar {
namespace multilevel {
static inline void partition(Hypergraph& hypergraph, const Context& context, const bool top_level);
}  // namespace multilevel
}  // namespace mt_kahypar

#include "mt-kahypar/partition/initial_partitioning/direct_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/recursive_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/recursive_bisection_initial_partitioner.h"

#define REGISTER_INITIAL_PARTITIONER(id, partitioner)                                  \
  static kahypar::meta::Registrar<InitialPartitionerFactory> register_ ## partitioner( \
    id,                                                                                \
    [](Hypergraph& hypergraph, const Context& context, const bool top_level)           \
    -> IInitialPartitioner* {                                                          \
    return new partitioner(hypergraph, context, top_level);                            \
  })

namespace mt_kahypar {
REGISTER_INITIAL_PARTITIONER(InitialPartitioningMode::direct, DirectInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningMode::recursive, RecursiveInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(InitialPartitioningMode::recursive_bisection, RecursiveBisectionInitialPartitioner);
}  // namespace mt_kahypar

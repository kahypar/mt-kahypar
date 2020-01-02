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
#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/partition/initial_partitioning/flat/random_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/flat/bfs_initial_partitioner.h"

#define REGISTER_FLAT_INITIAL_PARTITIONER(id, partitioner)                                           \
  static kahypar::meta::Registrar<FlatInitialPartitionerFactory> register_ ## partitioner(           \
    id,                                                                                              \
    [](tbb::task* parent, InitialPartitioningDataContainer& ip_hypergraph, const Context& context)   \
    -> tbb::task* {                                                                                  \
    return new(parent->allocate_child()) partitioner(ip_hypergraph, context);                        \
  })

namespace mt_kahypar {
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::random, RandomInitialPartitioner);
REGISTER_FLAT_INITIAL_PARTITIONER(InitialPartitioningAlgorithm::bfs, BFSInitialPartitioner);
}  // namespace mt_kahypar

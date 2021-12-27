/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/partition/initial_partitioning/deep_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/recursive_bipartitioning_initial_partitioner.h"

#define REGISTER_INITIAL_PARTITIONER(id, partitioner)                                           \
  static kahypar::meta::Registrar<InitialPartitionerFactory> register_ ## partitioner(          \
    id,                                                                                         \
    [](PartitionedHypergraph& hypergraph, const Context& context)                               \
    -> IInitialPartitioner* {                                                                   \
    return new partitioner(hypergraph, context);                                                \
  })

namespace mt_kahypar {
REGISTER_INITIAL_PARTITIONER(Mode::deep_multilevel, DeepInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(Mode::recursive_bipartitioning, RecursiveBipartitioningInitialPartitioner);
}  // namespace mt_kahypar

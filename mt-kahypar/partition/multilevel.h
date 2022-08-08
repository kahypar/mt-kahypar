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

#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar::multilevel {

// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB blocking-style.
PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context);
// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB continuation-style.
// ! Note, the final partitioned hypergraph is moved into the
// ! passed partitioned hypergraph object.
void partition_async(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context, tbb::task* parent);

// ! Performs a multilevel partitioning v-cycle on the given hypergraph
// ! in TBB blocking-style.
void partitionVCycle(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context);

}  // namespace mt_kahypar

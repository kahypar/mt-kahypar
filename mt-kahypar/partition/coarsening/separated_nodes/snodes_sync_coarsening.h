/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::SepNodesStack;
using ds::Array;

void coarsenSynchronized(SepNodesStack& stack, const Hypergraph& original_hg, const vec<Level>& levels,
                         const Context& context, const HypernodeID& start_num_nodes, const HypernodeID& target_num_nodes,
                         Array<PartitionID>* part_ids = nullptr);
} // namepace star_partitioning
} // namespace mt_kahypar

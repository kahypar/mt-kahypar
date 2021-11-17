/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/partition/context.h"


namespace mt_kahypar {
namespace recursive_bisection {

/*!
 * RECURSIVE BISECTION
 * The recursive bisection algorithm starts by performing a parallel multilevel bisection.
 * Once the hypergraph is bisected both blocks are partitioned recursively in parallel until the
 * desired number of blocks is reached.
*/

PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context);

void partition(PartitionedHypergraph& hypergraph, const Context& context);

}  // namespace recursive_bisection
}  // namespace mt_kahypar

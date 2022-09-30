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

#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/datastructures/static_graph_factory.h"
#include "mt-kahypar/datastructures/separated_nodes.h"

namespace mt_kahypar {

// forward
namespace ds {
template <typename Hypergraph, typename HypergraphFactory>
class PartitionedGraph;
}

class Context;
using PartitionedHypergraph = ds::PartitionedGraph<ds::StaticGraph, ds::StaticGraphFactory>;

namespace star_partitioning {

HyperedgeWeight rateSeparated(PartitionedHypergraph& phg, const Context& context, const HypernodeID u, const PartitionID to);

} // namepace star_partitioning
} // namepace mt_kahypar

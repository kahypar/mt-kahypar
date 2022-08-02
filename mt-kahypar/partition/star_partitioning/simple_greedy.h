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

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;

/*!
 * Simple star partitioning implementation that sorts all leaf nodes
 * and assigns them in a global greedy order.
*/
class SimpleGreedy {
 public:
  SimpleGreedy(const PartitionID& k): _k(k), _tmp_edge_weights(k) { }

  void partition(PartitionedHypergraph& hypergraph, const Context& context,
                 Array<HypernodeWeight>& part_weights, bool parallel = true);

 private:
  PartitionID _k;
  tbb::enumerable_thread_specific<Array<HyperedgeWeight>> _tmp_edge_weights;
};

} // namepace star_partitioning
} // namepace mt_kahypar


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
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;
using ds::SeparatedNodes;

HyperedgeWeight partition(PartitionedHypergraph& hypergraph, SeparatedNodes& s_nodes,
                          const Context& context, bool parallel = true);

void getEdgeWeightsOfNode(PartitionedHypergraph& phg, const SeparatedNodes& s_nodes, Array<HyperedgeWeight>& weights_per_part,
                          const HypernodeID& node, const HypernodeID& index = 0);

template<typename F, typename G>
auto compare_gain_weight_ratio(F get_gain, G get_node_weight) {
  return [get_gain, get_node_weight](const HypernodeID& left, const HypernodeID& right) {
    const HypernodeWeight weight_left = get_node_weight(left);
    const HypernodeWeight weight_right = get_node_weight(right);
    if (weight_left == 0) {
    return false;
    } else if (weight_right == 0) {
    return true;
    }
    return static_cast<double>(get_gain(left)) / weight_left
            < static_cast<double>(get_gain(right)) / weight_right;
  };
}

} // namepace star_partitioning
} // namepace mt_kahypar

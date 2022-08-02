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

#include <vector>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

#include <tbb/parallel_sort.h>

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;

class Approximate {
 public:
  Approximate(const PartitionID& k): _k(k) { }

  void partition(PartitionedHypergraph& phg, const Context& context,
                 Array<HypernodeWeight>& part_weights, parallel_tag_t);

  void partition(PartitionedHypergraph& phg, const Context& context,
                 Array<HypernodeWeight>& part_weights);

 private:
  PartitionID _k;
};

/*
* Returns the nodes that are _not_ included in ascending order.
*/
template<typename F, typename G>
std::vector<HypernodeID> minKnapsack(const parallel::scalable_vector<HypernodeID>& sorted_nodes,
                                     const HypernodeWeight& capacity,
                                     F get_gain_fn, G get_node_weight_fn) {
  HypernodeWeight total_weight = 0;
  for (const HypernodeID& node: sorted_nodes) {
    total_weight += get_node_weight_fn(node);
  }
  std::vector<HypernodeID> result;
  const HypernodeWeight excluded_weight = total_weight - capacity;
  if (excluded_weight <= 0) {
    return result;
  }

  HyperedgeWeight min_gain = std::numeric_limits<HyperedgeWeight>::max();
  HypernodeID min_element = kInvalidHypernode;
  size_t min_num_elements = 0;
  HypernodeWeight current_weight = 0;
  HyperedgeWeight current_gain = 0;
  for (const HypernodeID& node: sorted_nodes) {
    const HypernodeWeight weight = get_node_weight_fn(node);
    const double gain = static_cast<double>(get_gain_fn(node));
    if (current_weight + weight <= excluded_weight) {
      current_weight += weight;
      current_gain += gain;
      result.push_back(node);
      if (current_weight == excluded_weight && current_gain <= min_gain) {
        return result;
      }
    } else if (weight == 0 || current_gain + (gain / weight) * (excluded_weight - current_weight)
               > static_cast<double>(min_gain)) {
      // bounding
      break;
    } else if (current_gain + gain < min_gain) {
      min_gain = current_gain + gain;
      min_element = node;
      min_num_elements = result.size();
    }
  }

  // reconstruct result
  result.resize(min_num_elements);
  result.push_back(min_element);
  return result;
}

} // namepace star_partitioning
} // namepace mt_kahypar


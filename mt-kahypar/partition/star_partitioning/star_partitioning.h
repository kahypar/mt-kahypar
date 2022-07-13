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
#include "mt-kahypar/partition/star_partitioning/approximate.h"
#include "mt-kahypar/partition/star_partitioning/simple_greedy.h"

namespace mt_kahypar {
namespace star_partitioning {
using ds::Array;

template<typename F, typename G, typename H>
void partition(PartitionedHypergraph& hypergraph, const Context& context,
               const HypernodeID num_nodes, const std::vector<HypernodeWeight>& max_part_weights,
               F get_edge_weights_of_node_fn, G get_node_weight_fn, H set_part_id_fn) {
    Array<HypernodeWeight> part_weights(hypergraph.k());
    for (PartitionID part = 0; part < hypergraph.k(); ++part) {
        part_weights[part] = hypergraph.partWeight(part);
    }

    if (context.partition.star_partitioning_algorithm == StarPartitioningAlgorithm::simple_greedy) {
        SimpleGreedy sg(context.partition.k);
        sg.partition(num_nodes, part_weights, max_part_weights, get_edge_weights_of_node_fn,
                    get_node_weight_fn, set_part_id_fn);
    } else if (context.partition.star_partitioning_algorithm == StarPartitioningAlgorithm::approximate) {
        Approximate ap(context.partition.k);
        ap.partition(num_nodes, part_weights, max_part_weights, get_edge_weights_of_node_fn,
                    get_node_weight_fn, set_part_id_fn);
    }
}

} // namepace star_partitioning
} // namepace mt_kahypar

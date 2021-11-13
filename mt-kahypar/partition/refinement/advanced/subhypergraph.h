/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "definitions.h" // whfc::Node
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/shared_mapping.h"

namespace mt_kahypar {

using SharedMap = ds::SharedMapping<SearchID, HypernodeID, whfc::Node>;

struct Subhypergraph {
  PartitionID block_0;
  PartitionID block_1;
  vec<HypernodeID> nodes_of_block_0;
  vec<HypernodeID> nodes_of_block_1;
  HypernodeWeight weight_of_block_0;
  HypernodeWeight weight_of_block_1;
  vec<HyperedgeID> hes;
  size_t num_pins;

  size_t numNodes() const {
    return nodes_of_block_0.size() + nodes_of_block_1.size();
  }
};

inline std::ostream& operator<<(std::ostream& out, const Subhypergraph& sub_hg) {
  out << "[Nodes=" << sub_hg.numNodes()
      << ", Edges=" << sub_hg.hes.size()
      << ", Pins=" << sub_hg.num_pins
      << ", Blocks=(" << sub_hg.block_0 << "," << sub_hg.block_1 << ")"
      << ", Weights=(" << sub_hg.weight_of_block_0 << "," << sub_hg.weight_of_block_1 << ")]";
  return out;
}

}  // namespace kahypar

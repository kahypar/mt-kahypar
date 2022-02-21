/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "datastructure/flow_hypergraph_builder.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

enum class MoveSequenceState : uint8_t {
  IN_PROGRESS = 0,
  SUCCESS = 1,
  VIOLATES_BALANCE_CONSTRAINT = 2,
  WORSEN_SOLUTION_QUALITY = 3,
  WORSEN_SOLUTION_QUALITY_WITHOUT_ROLLBACK = 4,
  TIME_LIMIT = 5
};

// Represents a sequence of vertex moves with an
// expected improvement of the solution quality if we
// apply the moves
struct MoveSequence {
  vec<Move> moves;
  Gain expected_improvement; // >= 0
  MoveSequenceState state = MoveSequenceState::IN_PROGRESS;
};

struct FlowProblem {
  whfc::Node source;
  whfc::Node sink;
  HyperedgeWeight total_cut;
  HyperedgeWeight non_removable_cut;
  HypernodeWeight weight_of_block_0;
  HypernodeWeight weight_of_block_1;
};

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

} // namespace mt_kahypar
/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <cstdint>

#include "mt-kahypar/macros.h"

namespace mt_kahypar {

using RatingType = double;
using HypernodeID = uint64_t;
using HyperedgeID = uint64_t;
using HypernodeWeight = int32_t;
using HyperedgeWeight = int32_t;
using PartitionID = int32_t;
using Gain = HyperedgeWeight;

struct Move {
  PartitionID from;
  PartitionID to;
  Gain gain;
};

// Constant Declarations
static constexpr PartitionID kInvalidPartition = -1;
static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
static constexpr HypernodeID kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();
static constexpr size_t kEdgeHashSeed = 42;

namespace common {

static constexpr size_t NUMA_NODE_IDENTIFIER = 48;

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_global_vertex_id(const int node, const size_t vertex_pos) {
  return ( static_cast<HypernodeID>(node) << NUMA_NODE_IDENTIFIER ) | vertex_pos;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_position_of_vertex(const HypernodeID u) {
  return ((1UL << NUMA_NODE_IDENTIFIER) - 1) & u;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_vertex(const HypernodeID u) {
  return static_cast<int>(u >> NUMA_NODE_IDENTIFIER);
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_global_edge_id(const int node, const size_t edge_pos) {
  return ( static_cast<HypernodeID>(node) << NUMA_NODE_IDENTIFIER ) | edge_pos;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_local_position_of_edge(const HyperedgeID e) {
  return ((1UL << NUMA_NODE_IDENTIFIER) - 1) & e;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_edge(const HyperedgeID e) {
  return static_cast<int>(e >> NUMA_NODE_IDENTIFIER);
}

} // namespace common

} // namespace mt_kahypar

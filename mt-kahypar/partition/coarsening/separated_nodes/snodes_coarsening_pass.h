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

#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace star_partitioning {

using ds::SeparatedNodes;
using ds::Array;

enum class SNodesCoarseningStage : uint8_t {
  DEGREE_ZERO = 0,
  PREFERABLE_DEGREE_ONE = 1,
  DEGREE_ONE_AND_TWINS = 2,
  DEGREE_ONE_AND_DEGREE_TWO_AND_TWINS = 3,
  // TODO: mixing with high degree??
  ANY_DEGREE = 4,
  ANY_DEGREE_RELAXED = 5,
  ANYTHING = 6
};

bool allowsDegreeTwo(const SNodesCoarseningStage& stage);

bool allowsHighDegree(const SNodesCoarseningStage& stage);

bool appliesTwins(const SNodesCoarseningStage& stage);

class SNodesCoarseningPass {
  struct NodeInfo {
    NodeInfo(): density(0.0), assigned_graph_node(kInvalidHypernode) { }

    double density;
    HypernodeID assigned_graph_node;
    // TODO: position in spanning tree
  };

  struct FullNodeInfo {
    FullNodeInfo():
                 node(kInvalidHypernode),
                 degree(0),
                 density(0),
                 assigned_graph_node(kInvalidHypernode)  { }

    FullNodeInfo(HypernodeID node, HyperedgeID degree, const NodeInfo& info):
                 node(node),
                 degree(degree),
                 density(info.density),
                 assigned_graph_node(info.assigned_graph_node)  { }

    HypernodeID node;
    HyperedgeID degree;
    double density;
    HypernodeID assigned_graph_node;
    // TODO: position in spanning tree
  };

  struct LocalizedData {
    tbb::enumerable_thread_specific<HypernodeID> match_counter;
    tbb::enumerable_thread_specific<vec<HypernodeID>> degree_one_nodes;
    tbb::enumerable_thread_specific<vec<HypernodeID>> degree_two_nodes;
  };

  struct Params {
    double accepted_density_diff;
    HypernodeWeight max_node_weight;
    HypernodeID degree_one_cluster_size;
  };

  // tuning constants
  static const HypernodeID MAX_CLUSTER_SIZE = 4;
  static constexpr double PREFERRED_DENSITY_DIFF = 1.6;
  static constexpr double TOLERABLE_DENSITY_DIFF = 2.1;
  static constexpr double RELAXED_DENSITY_DIFF = 4.1;

 public:
  SNodesCoarseningPass(const Hypergraph& hg, const Context& context,
                       const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage);

  // note: communities are allowed to be an empty vec
  void run(vec<HypernodeID>& communities);

  // only for testing
  const Array<HypernodeID>& nodeInfoBegin() const {
    return _node_info_begin;
  }

  // only for testing
  const Array<FullNodeInfo>& nodeInfo() const {
    return _node_info;
  }

 private:
  void setupNodeInfo();

  HypernodeID runCurrentStage(vec<HypernodeID>& communities);

  void applyCoarseningForNode(const Params& params, vec<HypernodeID>& communities,
                              LocalizedData& data, const HypernodeID& node);

  void sortByDensity(vec<HypernodeID>& nodes);

  const FullNodeInfo& info(HypernodeID index) const {
    return _node_info[index];
  }

  const Hypergraph& _hg;
  const Context& _context;
  const SeparatedNodes& _s_nodes;
  Array<HypernodeID> _node_info_begin;
  Array<FullNodeInfo> _node_info;
  HypernodeID _current_num_nodes;
  HypernodeID _target_num_nodes;
  SNodesCoarseningStage _stage;
};

} // namepace star_partitioning
} // namespace mt_kahypar

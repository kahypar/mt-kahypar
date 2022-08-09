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
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace star_partitioning {

using ds::Array;

/*
 * Spanning tree that is built on the coarsened graph to provide a distance metric.
*/
class SpanningTree {
  static constexpr uint8_t kInvalidDepth = std::numeric_limits<uint8_t>::max();

 public:
  static const size_t max_children = 7;

  SpanningTree(HypernodeID num_nodes, HypernodeID root, uint8_t max_depth):
          _max_depth(max_depth), _num_nodes(num_nodes), _nodes(), _depth() {
    ASSERT(num_nodes < 50000);
    ASSERT(root < num_nodes);
    _nodes.assign(num_nodes * max_children, kInvalidHypernode);
    _depth.assign(num_nodes, kInvalidDepth);
    _depth[root] = 0;
  }

  bool isLeaf(HypernodeID node) {
    ASSERT(node < _depth.size() && depth(node) != kInvalidDepth);
    return _nodes[node * max_children] == kInvalidHypernode;
  }

  uint8_t depth(HypernodeID node) {
    ASSERT(node < _depth.size());
    return _depth[node];
  }

  uint8_t requiredDepthForChild(HypernodeID node) {
    auto [target_node, slot, do_allocate]  = locateSlotForChild(node);
    return depth(target_node) + (do_allocate ? 2 : 1);
  }

  bool mayAddChild(HypernodeID node) {
    ASSERT(node < _num_nodes);
    return requiredDepthForChild(node) < _max_depth;
  }

  void addChild(HypernodeID node, HypernodeID child);

  IteratorRange<const HypernodeID*> children(HypernodeID node) const {
    ASSERT(node * max_children < _nodes.size());
    size_t i = node * max_children;
    for (; i < (node + 1) * max_children; ++i) {
      if (_nodes[i] == kInvalidHypernode) {
        break;
      }
    }
    return IteratorRange<const HypernodeID*>(
      _nodes.data() + node * max_children,
      _nodes.data() + i);
  }

 private:
  // the node, the slot index of the node, whether a pseudo node needs to be allocated
  std::tuple<HypernodeID, HypernodeID, bool> locateSlotForChild(HypernodeID node);

  bool isPseudoNode(HypernodeID node) {
    ASSERT(node != kInvalidHypernode);
    return node >= _num_nodes;
  }

  uint8_t _max_depth;
  HypernodeID _num_nodes;
  vec<HypernodeID> _nodes;
  vec<uint8_t> _depth;
};

} // namespace mt_kahypar
} // namespace star_partitioning

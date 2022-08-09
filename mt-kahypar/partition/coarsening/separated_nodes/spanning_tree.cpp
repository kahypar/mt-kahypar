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

#include "spanning_tree.h"

namespace mt_kahypar {
namespace star_partitioning {

void SpanningTree::addChild(HypernodeID node, HypernodeID child) {
  ASSERT(node < _num_nodes && child < _num_nodes);
  auto [target_node, slot, do_allocate] = locateSlotForChild(node);
  if (do_allocate) {
    ASSERT(isLeaf(_nodes[target_node * max_children + slot]));
    const size_t first_slot_of_pseudo_node = _nodes.size();
    const HypernodeID pseudo_node = first_slot_of_pseudo_node / max_children;
    _nodes.resize(_nodes.size() + max_children, kInvalidHypernode);
    _depth.push_back(depth(target_node) + 1);

    const HypernodeID swapped_node = _nodes[target_node * max_children + slot];
    _nodes[target_node * max_children + slot] = pseudo_node;
    _nodes[first_slot_of_pseudo_node] = swapped_node;
    _nodes[first_slot_of_pseudo_node + 1] = child;
    _depth[pseudo_node] = depth(target_node) + 1;
    _depth[swapped_node] = depth(target_node) + 2;
    _depth[child] = depth(target_node) + 2;
  } else {
    ASSERT(_nodes[target_node * max_children + slot] == kInvalidHypernode);
    _nodes[target_node * max_children + slot] = child;
    _depth[child] = depth(target_node) + 1;
  }
  ASSERT(_depth[child] < _max_depth);
}

std::tuple<HypernodeID, HypernodeID, bool> SpanningTree::locateSlotForChild(HypernodeID node) {
  ASSERT(node < _depth.size());
  const uint8_t parent_depth = depth(node);
  const size_t start = node * max_children;
  const size_t end = (node + 1) * max_children;

  bool allocate_pseudo_node = false;
  std::tuple<HypernodeID, HypernodeID, bool> best;
  uint8_t best_depth = kInvalidDepth;
  for (size_t i = 0; start + i < end; ++i) {
    const HypernodeID current_node = _nodes[start + i];
    if (current_node == kInvalidHypernode) {
      return {node, i, false};
    } else if (isPseudoNode(current_node)) {
      ASSERT(depth(current_node) == parent_depth + 1);
      // if there is already an allocated pseudo node, will need to use one of those
      allocate_pseudo_node = true;
      auto [target_node, slot, do_allocate]  = locateSlotForChild(current_node);
      const uint8_t new_depth = depth(target_node) + do_allocate ? 2 : 1;
      if (new_depth == parent_depth + 1) {
        return {target_node, slot, do_allocate};
      } else if (new_depth < best_depth) {
        best = {target_node, slot, do_allocate};
        best_depth = new_depth;
      }
    } else if (allocate_pseudo_node && isLeaf(current_node)) {
      ASSERT(depth(current_node) == parent_depth + 1);
      return {node, i, true};
    }
  }
  ASSERT(false);
  return {0, 0, false};
}

} // namespace mt_kahypar
} // namespace star_partitioning

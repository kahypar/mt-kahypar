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

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/datastructures/priority_queue.h"

namespace mt_kahypar {
namespace star_partitioning {

void SpanningTree::addChild(HypernodeID node, HypernodeID child) {
  ASSERT(child < _num_nodes);
  ASSERT(contains(node) && !contains(child));
  auto [target_node, slot, do_allocate] = locateSlotForChild(node);
  if (do_allocate) {
    ASSERT(_nodes[target_node * max_children + slot] == kInvalidHypernode
           || isLeaf(_nodes[target_node * max_children + slot]));
    const size_t first_slot_of_pseudo_node = _nodes.size();
    const HypernodeID pseudo_node = first_slot_of_pseudo_node / max_children;
    _nodes.resize(_nodes.size() + max_children, kInvalidHypernode);
    _depth.push_back(depth(target_node) + 1);

    const HypernodeID swapped_node = _nodes[target_node * max_children + slot];
    _nodes[target_node * max_children + slot] = pseudo_node;
    _nodes[first_slot_of_pseudo_node] = child;
    _nodes[first_slot_of_pseudo_node + 1] = swapped_node;
    _depth[pseudo_node] = depth(target_node) + 1;
    _depth[child] = depth(target_node) + 2;
    if (swapped_node != kInvalidHypernode) {
      _depth[swapped_node] = depth(target_node) + 2;
    }
  } else {
    ASSERT(_nodes[target_node * max_children + slot] == kInvalidHypernode);
    _nodes[target_node * max_children + slot] = child;
    _depth[child] = depth(target_node) + 1;
  }
  ASSERT(_depth[child] < _max_depth);
}

std::tuple<HypernodeID, HypernodeID, bool> SpanningTree::locateSlotForChild(HypernodeID node) const {
  ASSERT(contains(node));
  const uint8_t parent_depth = depth(node);
  const size_t start = node * max_children;
  const size_t end = (node + 1) * max_children - 1;
  const HypernodeID last_node = _nodes[end];

  const bool allocate_pseudo_node = (last_node != kInvalidHypernode);
  if (allocate_pseudo_node) {
    ASSERT(isPseudoNode(last_node));
    std::tuple<HypernodeID, HypernodeID, bool> best;
    uint8_t best_depth = kInvalidDepth;
    for (size_t i = 0; i < max_children; ++i) {
      const HypernodeID current_node = _nodes[end - i];
      ASSERT(current_node != kInvalidHypernode);
      if (isPseudoNode(current_node)) {
        // if there is already an allocated pseudo node, will need to use one of those
        auto [target_node, slot, do_allocate] = locateSlotForChild(current_node);
        const uint8_t new_depth = depth(target_node) + (do_allocate ? 2 : 1);
        if (new_depth == parent_depth + 1) {
          return {target_node, slot, do_allocate};
        } else if (new_depth < best_depth) {
          best = {target_node, slot, do_allocate};
          best_depth = new_depth;
        }
      } else if (isLeaf(current_node) && (best_depth > parent_depth + 2)) {
        return {node, max_children - i - 1, true};
      }
    }
    return best;
  } else {
    for (size_t i = 0; start + i < end; ++i) {
      if (_nodes[start + i] == kInvalidHypernode) {
        return {node, i, false};
      }
    }
    ASSERT(last_node == kInvalidHypernode);
    return {node, max_children - 1, true};
  }
}

SpanningTree constructMaxSpanningTree(const Hypergraph& hg, const HypernodeID& max_depth) {
  using PriorityQueue = ds::ExclusiveHandleHeap< ds::MaxHeap<Gain, PartitionID> >;

  Array<HypernodeID> node_to_current_parent;
  PriorityQueue pq(0);
  HypernodeID root;
  SpanningTree tree;
  tbb::parallel_invoke([&] {
    node_to_current_parent.assign(hg.initialNumNodes(), kInvalidHypernode);
  }, [&] {
    pq = std::move(PriorityQueue(hg.initialNumNodes()));
  }, [&] {
    tbb::enumerable_thread_specific<std::pair<HypernodeID, HypernodeWeight>>
              high_degree_node(std::make_pair(kInvalidHypernode, -1));
    tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](const HypernodeID& node) {
      auto& [best_node, best_incident_weight] = high_degree_node.local();
      if (hg.incidentWeight(node) > best_incident_weight) {
        best_incident_weight = hg.incidentWeight(node);
        best_node = node;
      }
    });
    root = high_degree_node.combine([&](auto left, auto right) {
      return left.second > right.second ? left : right;
    }).first;
    tree = SpanningTree(hg.initialNumNodes(), root, max_depth);
  });

  // use modified jarnik prim algorithm
  for (const HyperedgeID e: hg.incidentEdges(root)) {
    pq.insert(hg.edgeTarget(e), hg.edgeWeight(e));
    node_to_current_parent[hg.edgeTarget(e)] = root;
  }
  while (!pq.empty()) {
    const HypernodeID current = pq.top();
    const HypernodeID parent = node_to_current_parent[current];
    ASSERT(!tree.contains(current) && parent != kInvalidHypernode);
    pq.deleteTop();
    if (tree.mayAddChild(parent)) {
      tree.addChild(parent, current);
      if (tree.mayAddChild(current)) {
        for (const HyperedgeID e: hg.incidentEdges(current)) {
          const HypernodeID next = hg.edgeTarget(e);
          if (!tree.contains(next) && !pq.contains(next)) {
            pq.insert(next, hg.edgeWeight(e));
            node_to_current_parent[next] = current;
          } else if (!tree.contains(next) && hg.edgeWeight(e) > pq.getKey(next)) {
            pq.increaseKey(next, hg.edgeWeight(e));
            node_to_current_parent[next] = current;
          }
        }
      }
    } else {
      // reinsert with valid edge, if possible
      // (otherwise the node won't be part of the tree)
      HypernodeID best_parent = kInvalidHypernode;
      HypernodeWeight best_weight = -1;
      for (const HyperedgeID e: hg.incidentEdges(current)) {
        const HypernodeID possible_parent = hg.edgeTarget(e);
        if (possible_parent != parent
            && tree.contains(possible_parent)
            && tree.mayAddChild(possible_parent)
            && hg.edgeWeight(e) > best_weight) {
          best_parent = possible_parent;
          best_weight = hg.edgeWeight(e);
        }
      }
      if (best_parent != kInvalidHypernode) {
        pq.insert(current, best_weight);
        node_to_current_parent[current] = best_parent;
      }
    }
  }
  return tree;
}

} // namespace mt_kahypar
} // namespace star_partitioning

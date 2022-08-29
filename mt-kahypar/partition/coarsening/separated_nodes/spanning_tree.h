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
#include "mt-kahypar/utils/bit_ops.h"

namespace mt_kahypar {
namespace star_partitioning {

using ds::Array;

/*
 * Represents a path in a (depth-constrained) tree using a 64 bit integer.
*/
class TreePath {
  using ValT = uint64_t;
  static_assert(sizeof(ValT) * 8 == 64);

  static constexpr uint32_t BITS_PER_STEP = 3;
  static constexpr uint32_t MAX_STEP = (1UL << BITS_PER_STEP);
  static constexpr ValT LEADING_ONE = (1UL << 63);

 public:
  static constexpr uint32_t max_path_length = 63 / BITS_PER_STEP;

  explicit TreePath(): _val(LEADING_ONE) { }

  explicit TreePath(ValT val): _val(val) {
    ASSERT(_val == 0 || (_val | LEADING_ONE) != 0);
  }

  bool isValid() const {
    ASSERT(_val == 0 || (_val | LEADING_ONE) != 0);
    return _val != 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE uint32_t length() const {
    ASSERT(isValid());
    int used_bits = 63 - utils::lowest_set_bit_64(_val); // highest bit is reserved
    return (used_bits + BITS_PER_STEP - 1) / BITS_PER_STEP;
  }

  void append(uint32_t step) {
    ASSERT(step + 1 < MAX_STEP);
    ASSERT(isValid() && length() + 1 <= max_path_length);
    uint32_t n_zeros = 63 - BITS_PER_STEP;
    uint32_t shift = n_zeros - length() * BITS_PER_STEP;
    // zeros indicate end of path, thus add one to step
    _val |= (static_cast<ValT>(step + 1) << shift);
  }

  uint32_t distance(const TreePath& other) const {
    if (!isValid() || !other.isValid()) {
      return std::numeric_limits<uint32_t>::max();
    }

    ValT diff = _val ^ other._val;
    if (diff == 0) {
      return 0;
    }
    uint32_t common_length = (utils::leading_zeros_64(diff) - 1) / BITS_PER_STEP;
    ASSERT(common_length <= length() && common_length <= other.length());
    return length() + other.length() - 2 * common_length;
  }

  // comparison operators
  bool operator==(const TreePath& other) const {
    return _val == other._val;
  }

  bool operator!=(const TreePath& other) const {
    return _val != other._val;
  }

  bool operator<(const TreePath& other) const {
    return _val < other._val;
  }

  bool operator>(const TreePath& other) const {
    return _val > other._val;
  }

  bool operator<=(const TreePath& other) const {
    return _val <= other._val;
  }

  bool operator>=(const TreePath& other) const {
    return _val >= other._val;
  }

 private:
  ValT _val;
};

static_assert(std::is_trivially_copyable<TreePath>::value);
static const TreePath kInvalidPath = TreePath(0);


/*
 * Spanning tree that is built on the coarsened graph to provide a distance metric.
*/
class SpanningTree {
  static constexpr uint8_t kInvalidDepth = std::numeric_limits<uint8_t>::max();

 public:
  static const size_t max_children = 7;

  explicit SpanningTree(): _max_depth(0), _num_nodes(0), _root(kInvalidHypernode), _nodes(), _depth() { }

  explicit SpanningTree(HypernodeID num_nodes, HypernodeID root, uint8_t max_depth):
          _max_depth(max_depth), _num_nodes(num_nodes), _root(root), _nodes(), _depth() {
    ASSERT(num_nodes < 50000);
    ASSERT(root < num_nodes);
    _nodes.assign(num_nodes * max_children, kInvalidHypernode);
    _depth.assign(num_nodes, kInvalidDepth);
    _depth[root] = 0;
  }

  bool contains(HypernodeID node) const {
    ASSERT(node < _depth.size());
    return _depth[node] != kInvalidDepth;
  }

  bool isLeaf(HypernodeID node) const {
    ASSERT(node < _depth.size() && contains(node));
    return _nodes[node * max_children] == kInvalidHypernode;
  }

  uint8_t depth(HypernodeID node) const {
    ASSERT(node < _depth.size());
    return _depth[node];
  }

  uint8_t requiredDepthForChild(HypernodeID node) const {
    ASSERT(contains(node));
    auto [target_node, slot, do_allocate]  = locateSlotForChild(node);
    return depth(target_node) + (do_allocate ? 2 : 1);
  }

  bool mayAddChild(HypernodeID node) const {
    ASSERT(contains(node));
    return requiredDepthForChild(node) < _max_depth;
  }

  template<typename F>
  void calculatePaths(F receive_path) const {
    TreePath root_path;
    calculatePathsImpl(root_path, _root, receive_path);
  }

  void addChild(HypernodeID node, HypernodeID child);

  IteratorRange<const HypernodeID*> children(HypernodeID node) const {
    ASSERT(contains(node));
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
  std::tuple<HypernodeID, HypernodeID, bool> locateSlotForChild(HypernodeID node) const;

  bool isPseudoNode(HypernodeID node) const {
    ASSERT(node != kInvalidHypernode);
    return node >= _num_nodes;
  }

  template<typename F>
  void calculatePathsImpl(TreePath path, HypernodeID node, F receive_path) const {
    if (!isPseudoNode(node)) {
      receive_path(node, path);
    }
    uint32_t step = 0;
    for (const HypernodeID& child: children(node)) {
      TreePath child_path = path;
      child_path.append(step);
      calculatePathsImpl(child_path, child, receive_path);
      step++;
    }
  }

  uint8_t _max_depth;
  HypernodeID _num_nodes;
  HypernodeID _root;
  vec<HypernodeID> _nodes;
  vec<uint8_t> _depth;
};

SpanningTree constructMaxSpanningTree(const Hypergraph& hg, const HypernodeID& max_depth);

} // namespace mt_kahypar
} // namespace star_partitioning
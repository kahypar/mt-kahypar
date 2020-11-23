/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <tbb/parallel_for.h>

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/utils/range.h"

namespace mt_kahypar {
namespace ds {

class ContractionTree {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr size_t kInvalidVersion = std::numeric_limits<size_t>::max();

  using Timepoint = HypernodeID;

 public:
  struct Interval {
    explicit Interval() :
      start(kInvalidHypernode),
      end(kInvalidHypernode) { }

    Timepoint start;
    Timepoint end;
  };

 private:
  /**
   * Represents a node in contraction tree and contains all information
   * associated with that node.
   */
  class Node {
    public:
      Node() :
        _parent(0),
        _pending_contractions(0),
        _subtree_size(0),
        _version(kInvalidVersion),
        _interval() { }

      inline HypernodeID parent() const {
        return _parent;
      }

      inline void setParent(const HypernodeID parent) {
        _parent = parent;
      }

      inline HypernodeID pendingContractions() const {
        return _pending_contractions;
      }

      inline void incrementPendingContractions() {
        ++_pending_contractions;
      }

      inline void decrementPendingContractions() {
        --_pending_contractions;
      }

      inline HypernodeID subtreeSize() const {
        return _subtree_size;
      }

      inline void setSubtreeSize(const HypernodeID subtree_size) {
        _subtree_size = subtree_size;
      }

      inline size_t version() const {
        return _version;
      }

      inline void setVersion(const size_t version) {
        _version = version;
      }

      inline Interval interval() const {
        return _interval;
      }

      inline void setInterval(const Timepoint start, const Timepoint end) {
        ASSERT(start < end);
        _interval.start = start;
        _interval.end = end;
      }

      inline void reset(const HypernodeID u) {
        _parent = u;
        _pending_contractions = 0;
        _subtree_size = 0;
        _version = kInvalidVersion;
        _interval.start = kInvalidHypernode;
        _interval.end = kInvalidHypernode;
      }

    private:
      // ! Parent in the contraction tree
      HypernodeID _parent;
      // ! Number of pending contractions
      HypernodeID _pending_contractions;
      // ! Size of the subtree
      HypernodeID _subtree_size;
      // ! Version number of the hypergraph for which contract the corresponding vertex
      size_t _version;
      // ! "Time" interval on which the contraction of this node takes place
      Interval _interval;
  };

  static_assert(std::is_trivially_copyable<Node>::value, "Node is not trivially copyable");

 public:
  // ! Iterator to iterate over the childs of a tree node
  using ChildIterator = typename parallel::scalable_vector<HypernodeID>::const_iterator;

  explicit ContractionTree() :
    _num_hypernodes(0),
    _finalized(false),
    _tree(),
    _roots(),
    _version_roots(),
    _out_degrees(),
    _incidence_array() { }

  ContractionTree(ContractionTree&& other) :
    _num_hypernodes(other._num_hypernodes),
    _finalized(other._finalized),
    _tree(std::move(other._tree)),
    _roots(std::move(other._roots)),
    _version_roots(std::move(other._version_roots)),
    _out_degrees(std::move(other._out_degrees)),
    _incidence_array(std::move(other._incidence_array)) { }

  ContractionTree& operator= (ContractionTree&& other) {
    _num_hypernodes = other._num_hypernodes;
    _finalized = other._finalized;
    _tree = std::move(other._tree);
    _roots = std::move(other._roots);
    _version_roots = std::move(other._version_roots);
    _out_degrees = std::move(other._out_degrees);
    _incidence_array = std::move(other._incidence_array);
    return *this;
  }

  ~ContractionTree() {
    freeInternalData();
  }

  // ####################### Tree Node Information #######################

  HypernodeID num_hypernodes() const {
    return _num_hypernodes;
  }

  // ! Returns the parent of node u
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID parent(const HypernodeID u) const {
    return node(u).parent();
  }

  // ! Number of pending contractions of node u
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID pendingContractions(const HypernodeID u) const {
    return node(u).pendingContractions();
  }

  // ! Subtree size of node u
  HypernodeID subtreeSize(const HypernodeID u) const {
    ASSERT(_finalized, "Information currently not available");
    return node(u).subtreeSize();
  }


  size_t version(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _tree[u].version();
  }

  // ! Degree/Number of childs of node u
  HypernodeID degree(const HypernodeID u) const {
    ASSERT(_finalized, "Information currently not available");
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _out_degrees[u + 1] - _out_degrees[u];
  }

  const parallel::scalable_vector<HypernodeID>& roots() const {
    return _roots;
  }

  const parallel::scalable_vector<HypernodeID>& roots_of_version(const size_t version) const {
    ASSERT(version < _version_roots.size());
    return _version_roots[version];
  }

  Interval interval(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return node(u).interval();
  }

  // ####################### Iterators #######################

  // ! Returns a range to loop over the childs of a tree node u.
  IteratorRange<ChildIterator> childs(const HypernodeID u) const {
    ASSERT(_finalized, "Information currently not available");
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<ChildIterator>(
      _incidence_array.cbegin() + _out_degrees[u],
      _incidence_array.cbegin() + _out_degrees[u + 1]);
  }

  // ! Calls function f for each child of vertex u with the corresponding version
  template<typename F>
  void doForEachChildOfVersion(const HypernodeID u, const size_t version, const F& f) const {
    ASSERT(_finalized, "Information currently not available");
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    for ( const HypernodeID& v : childs(u) ) {
      if ( _tree[v].version() == version ) {
        f(v);
      }
    }
  }

  // ####################### Contraction Functions #######################

  // ! Registers a contraction in the contraction tree
  void registerContraction(const HypernodeID u, const HypernodeID v, const size_t version = 0) {
    node(u).incrementPendingContractions();
    node(v).setParent(u);
    node(v).setVersion(version);
  }

  // ! Unregisters a contraction in the contraction tree
  void unregisterContraction(const HypernodeID u, const HypernodeID v,
                             const Timepoint start, const Timepoint end,
                             const bool failed = false) {
    ASSERT(node(v).parent() == u, "Node" << u << "is not parent of node" << v);
    ASSERT(node(u).pendingContractions() > 0, "There are no pending contractions for node" << u);
    node(u).decrementPendingContractions();
    if ( failed ) {
      node(v).setParent(v);
      node(v).setVersion(kInvalidVersion);
    } else {
      node(v).setInterval(start, end);
    }
  }

  // ! Only for testing
  void setParent(const HypernodeID u, const HypernodeID v, const size_t version = 0) {
    node(u).setParent(v);
    node(u).setVersion(version);
  }


  // ! Only for testing
  void decrementPendingContractions(const HypernodeID u) {
    node(u).decrementPendingContractions();
  }

  // ####################### Initialize / Finalize #######################

  // ! Initializes the data structure in parallel
  void initialize(const HypernodeID num_hypernodes);

  // ! Finalizes the contraction tree which involve reversing the parent pointers
  // ! such that the contraction tree can be traversed in a top-down fashion and
  // ! computing the subtree sizes.
  void finalize(const size_t num_versions = 1);

  // ####################### Copy #######################

  // ! Copy contraction tree in parallel
  ContractionTree copy(const TaskGroupID);

  // ! Copy contraction tree sequentially
  ContractionTree copy();

  // ! Resets internal data structures
  void reset();

  // ! Free internal data in parallel
  void freeInternalData();

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

 private:
  // ! Accessor for contraction tree-related information
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Node& node(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _tree[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Node& node(const HypernodeID u) {
    return const_cast<Node&>(static_cast<const ContractionTree&>(*this).node(u));
  }

  HypernodeID _num_hypernodes;
  bool _finalized;
  parallel::scalable_vector<Node> _tree;
  parallel::scalable_vector<HypernodeID> _roots;
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> _version_roots;
  parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeID>> _out_degrees;
  parallel::scalable_vector<HypernodeID> _incidence_array;
};

}  // namespace ds
}  // namespace mt_kahypar

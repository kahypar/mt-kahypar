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
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

class ContractionTree {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  /**
   * Represents a node in contraction tree and contains all information
   * associated with that node.
   */
  class Node {
    public:
      Node() :
        _parent(0),
        _pending_contractions(0),
        _subtree_size(0) { }

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

    private:
      // ! Parent in the contraction tree
      HypernodeID _parent;
      // ! Number of pending contractions
      HypernodeID _pending_contractions;
      // ! Size of the subtree
      HypernodeID _subtree_size;
  };

  static_assert(std::is_trivially_copyable<Node>::value, "Node is not trivially copyable");


  // ! Iterator to iterate over the childs of a tree node
  using ChildIterator = typename parallel::scalable_vector<HypernodeID>::const_iterator;

 public:
  explicit ContractionTree() :
    _num_hypernodes(0),
    _finalized(false),
    _tree(),
    _out_degrees(),
    _incidence_array() { }

  ContractionTree(ContractionTree&& other) :
    _num_hypernodes(other._num_hypernodes),
    _finalized(other._finalized),
    _tree(std::move(other._tree)),
    _out_degrees(std::move(other._out_degrees)),
    _incidence_array(std::move(other._incidence_array)) { }

  ContractionTree& operator= (ContractionTree&& other) {
    _num_hypernodes = other._num_hypernodes;
    _finalized = other._finalized;
    _tree = std::move(other._tree);
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
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID parent(const HypernodeID u) const {
    return node(u).parent();
  }

  // ! Number of pending contractions of node u
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeID pendingContractions(const HypernodeID u) const {
    return node(u).pendingContractions();
  }

  // ! Subtree size of node u
  HypernodeID subtreeSize(const HypernodeID u) const {
    ASSERT(_finalized, "Information currently not available");
    return node(u).subtreeSize();
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

  // ####################### Iterators #######################

  // ! Returns a range to loop over the childs of a tree node u.
  IteratorRange<ChildIterator> childs(const HypernodeID u) const {
    ASSERT(_finalized, "Information currently not available");
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<ChildIterator>(
      _incidence_array.cbegin() + _out_degrees[u],
      _incidence_array.cbegin() + _out_degrees[u + 1]);
  }

  // ####################### Contraction Functions #######################

  // ! Registers a contraction in the contraction tree
  void registerContraction(const HypernodeID u, const HypernodeID v) {
    node(u).incrementPendingContractions();
    node(v).setParent(u);
  }

  // ! Unregisters a contraction in the contraction tree
  void unregisterContraction(const HypernodeID u, const HypernodeID v, const bool failed = false) {
    node(u).decrementPendingContractions();
    if ( failed ) {
      node(v).setParent(v);
    }
  }

  // ! Only for testing
  void setParent(const HypernodeID u, const HypernodeID v) {
    node(u).setParent(v);
  }


  // ! Only for testing
  void decrementPendingContractions(const HypernodeID u) {
    node(u).decrementPendingContractions();
  }

  // ####################### Initialize / Finalize #######################

  // ! Initializes the data structure in parallel
  void initialize(const HypernodeID num_hypernodes) {
    _num_hypernodes = num_hypernodes;
    tbb::parallel_invoke([&] {
      _tree.resize(_num_hypernodes);
      tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID hn) {
        node(hn).setParent(hn);
      });
    }, [&] {
      _out_degrees.assign(_num_hypernodes + 1, parallel::IntegralAtomicWrapper<HypernodeID>(0));
    }, [&] {
      _incidence_array.resize(_num_hypernodes);
    });
  }

  // ! Finalizes the contraction tree which involve reversing the parent pointers
  // ! such that the contraction tree can be traversed in a top-down fashion and
  // ! computing the subtree sizes.
  void finalize() {
    ASSERT(!_finalized, "Contraction tree already finalized");
    // Compute out degrees of each tree node
    utils::Timer::instance().start_timer("compute_out_degrees", "Compute Out-Degree of Tree Nodes");
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID hn) {
      ASSERT(node(hn).pendingContractions() == 0, "There are pending contractions for node" << hn);
      const HypernodeID parent = node(hn).parent();
      if ( parent != hn ) {
        ASSERT(parent + 1 <= _num_hypernodes, "Parent" << parent << "does not exist!");
        ++_out_degrees[parent + 1];
      }
    });
    utils::Timer::instance().stop_timer("compute_out_degrees");

    // Compute prefix sum over out degrees which will be the index pointer into the incidence array
    utils::Timer::instance().start_timer("out_degree_prefix_sum", "Compute Out-Degree Prefix Sum");
    parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeID>> incidence_array_pos;
    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<HypernodeID>, parallel::scalable_vector>
      out_degree_prefix_sum(_out_degrees);
    tbb::parallel_invoke([&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _out_degrees.size()), out_degree_prefix_sum);
    }, [&] {
      incidence_array_pos.assign(_num_hypernodes, parallel::IntegralAtomicWrapper<HypernodeID>(0));
    });
    utils::Timer::instance().stop_timer("out_degree_prefix_sum");

    // Reverse parent pointer of contraction tree such that it can be traversed in top-down fashion
    utils::Timer::instance().start_timer("reverse_contraction_tree", "Reverse Contraction Tree");
    StreamingVector<HypernodeID> tmp_roots;
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID hn) {
      const HypernodeID parent = node(hn).parent();
      if ( parent != hn ) {
        const HypernodeID pos = _out_degrees[parent] + incidence_array_pos[parent]++;
        ASSERT(pos < _out_degrees[parent + 1]);
        _incidence_array[pos] = hn;
      } else {
        // In that case node hn is a root
        const bool contains_subtree = (_out_degrees[hn + 1] - _out_degrees[hn]) > 0;
        if ( contains_subtree ) {
          tmp_roots.stream(hn);
        }
      }
    });
    _roots = tmp_roots.copy_parallel();
    utils::Timer::instance().stop_timer("reverse_contraction_tree");

    _finalized = true;

    // Compute subtree sizes of each root in parallel via dfs
    utils::Timer::instance().start_timer("compute_subtree_sizes", "Compute Subtree Sizes");
    tbb::parallel_for(0UL, _roots.size(), [&](const size_t i) {
      parallel::scalable_vector<HypernodeID> dfs;
      dfs.push_back(_roots[i]);
      while( !dfs.empty() ) {
        const HypernodeID u = dfs.back();
        if ( subtreeSize(u) == 0 ) {
          // Visit u for the first time => push all childs on the dfs stack
          for ( const HypernodeID& v : childs(u)) {
            dfs.push_back(v);
          }
          // Mark u as visited
          node(u).setSubtreeSize(1);
        } else {
          // Visit u for second time => accumulate subtree sizes and pop u
          dfs.pop_back();
          HypernodeID subtree_size = 0;
          for ( const HypernodeID& v : childs(u) ) {
            subtree_size += ( subtreeSize(v) + 1 );
          }
          node(u).setSubtreeSize(subtree_size);
        }
      }
    });
    utils::Timer::instance().stop_timer("compute_subtree_sizes");

    tbb::parallel_invoke([&] {
      parallel::free(incidence_array_pos);
    }, [&] {
      tmp_roots.clear_parallel();
    });
  }

  // ####################### Copy #######################

  // ! Copy contraction tree in parallel
  ContractionTree copy(const TaskGroupID) {
    ContractionTree tree;

    tree._num_hypernodes = _num_hypernodes;
    tree._finalized = _finalized;

    tbb::parallel_invoke([&] {
      tree._tree.resize(_tree.size());
      memcpy(tree._tree.data(), _tree.data(),
        sizeof(Node) * _tree.size());
    }, [&] {
      tree._roots.resize(_roots.size());
      memcpy(tree._roots.data(), _roots.data(),
        sizeof(HypernodeID) * _roots.size());
    }, [&] {
      tree._out_degrees.resize(_out_degrees.size());
      for ( size_t i = 0; i < _out_degrees.size(); ++i ) {
        tree._out_degrees[i] = _out_degrees[i];
      }
    }, [&] {
      tree._incidence_array.resize(_incidence_array.size());
      memcpy(tree._incidence_array.data(), _incidence_array.data(),
        sizeof(HypernodeID) * _incidence_array.size());
    });

    return tree;
  }

  // ! Copy contraction tree sequentially
  ContractionTree copy() {
    ContractionTree tree;

    tree._num_hypernodes = _num_hypernodes;
    tree._finalized = _finalized;

    tree._tree.resize(_tree.size());
    memcpy(tree._tree.data(), _tree.data(),
      sizeof(Node) * _tree.size());
    tree._roots.resize(_roots.size());
    memcpy(tree._roots.data(), _roots.data(),
      sizeof(HypernodeID) * _roots.size());
    tree._out_degrees.resize(_out_degrees.size());
    for ( size_t i = 0; i < _out_degrees.size(); ++i ) {
      tree._out_degrees[i] = _out_degrees[i];
    }
    tree._incidence_array.resize(_incidence_array.size());
    memcpy(tree._incidence_array.data(), _incidence_array.data(),
      sizeof(HypernodeID) * _incidence_array.size());

    return tree;
  }

  // ! Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 ) {
      parallel::parallel_free(_tree, _roots, _out_degrees, _incidence_array);
    }
    _num_hypernodes = 0;
    _finalized = false;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Tree Nodes", sizeof(Node) * _tree.size());
    parent->addChild("Roots", sizeof(HypernodeID) * _roots.size());
    parent->addChild("Out-Degrees", sizeof(HypernodeID) * _out_degrees.size());
    parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
  }

 private:
  // ! Accessor for contraction tree-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Node& node(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _tree[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Node& node(const HypernodeID u) {
    return const_cast<Node&>(static_cast<const ContractionTree&>(*this).node(u));
  }

  HypernodeID _num_hypernodes;
  bool _finalized;
  parallel::scalable_vector<Node> _tree;
  parallel::scalable_vector<HypernodeID> _roots;
  parallel::scalable_vector<parallel::IntegralAtomicWrapper<HypernodeID>> _out_degrees;
  parallel::scalable_vector<HypernodeID> _incidence_array;
};

}  // namespace ds
}  // namespace mt_kahypar

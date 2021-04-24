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

#include "mt-kahypar/datastructures/contraction_tree.h"

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {
namespace ds {

// ! Initializes the data structure in parallel
void ContractionTree::initialize(const HypernodeID num_hypernodes) {
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
void ContractionTree::finalize(const size_t num_versions) {
  ASSERT(!_finalized, "Contraction tree already finalized");
  // Compute out degrees of each tree node
  utils::Timer::instance().start_timer("compute_out_degrees", "Compute Out-Degree of Tree Nodes");
  tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID hn) {
    ASSERT(node(hn).pendingContractions() == 0, "There are"
      << node(hn).pendingContractions() << "pending contractions for node" << hn);
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

  // Compute roots for each version
  // Each contraction/edge in the contraction tree is associated with a version.
  // Later we want to be able to traverse the contraction tree for a specific version
  // in a top-down fashion. Therefore, we compute for each version the corresponding roots.
  // A vertex is a root of a version if contains a child with that version less than
  // the version number of the vertex itself. Note, that for all vertices in the contraction
  // tree version(u) <= version(parent(u)).
  utils::Timer::instance().start_timer("compute_version_roots", "Compute Roots for each Version");
  parallel::scalable_vector<StreamingVector<HypernodeID>> tmp_version_roots(num_versions);
  tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID u) {
    std::sort(_incidence_array.begin() + _out_degrees[u],
              _incidence_array.begin() + _out_degrees[u + 1],
              [&](const HypernodeID& u, const HypernodeID& v) {
                const size_t u_version = version(u);
                const size_t v_version = version(v);
                const Interval& u_ival = node(u).interval();
                const Interval& v_ival = node(v).interval();
                return u_version < v_version ||
                  ( u_version == v_version && u_ival.end > v_ival.end ) ||
                  ( u_version == v_version && u_ival.end == v_ival.end && u_ival.start > v_ival.start ) ||
                  ( u_version == v_version && u_ival.end == v_ival.end && u_ival.start == v_ival.start && u < v );
              });

    size_t version_u = _tree[u].version();
    ASSERT(version_u <= _tree[_tree[u].parent()].version());
    size_t last_version = kInvalidVersion;
    for ( const HypernodeID& v : childs(u) ) {
      size_t version_v = _tree[v].version();
      ASSERT(version_v < num_versions, V(version_v) << V(num_versions));
      if ( version_v != last_version && version_v < version_u ) {
        tmp_version_roots[version_v].stream(u);
      }
      last_version = version_v;
    }
  });
  _version_roots.resize(num_versions);
  tbb::parallel_for(0UL, num_versions, [&](const size_t i) {
    _version_roots[i] = tmp_version_roots[i].copy_parallel();
    tmp_version_roots[i].clear_parallel();
  });
  utils::Timer::instance().stop_timer("compute_version_roots");

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
ContractionTree ContractionTree::copy(const TaskGroupID) {
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
    const size_t num_versions = _version_roots.size();
    tree._version_roots.resize(num_versions);
    tbb::parallel_for(0UL, num_versions, [&](const size_t i) {
      tree._version_roots[i].resize(_version_roots[i].size());
      memcpy(tree._version_roots[i].data(), _version_roots[i].data(),
        sizeof(HypernodeID) * _version_roots[i].size());
    });
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
ContractionTree ContractionTree::copy() {
  ContractionTree tree;

  tree._num_hypernodes = _num_hypernodes;
  tree._finalized = _finalized;

  tree._tree.resize(_tree.size());
  memcpy(tree._tree.data(), _tree.data(),
    sizeof(Node) * _tree.size());
  tree._roots.resize(_roots.size());
  memcpy(tree._roots.data(), _roots.data(),
    sizeof(HypernodeID) * _roots.size());
  const size_t num_versions = _version_roots.size();
  tree._version_roots.resize(num_versions);
  for ( size_t i = 0; i < num_versions; ++i ) {
    tree._version_roots[i].resize(_version_roots[i].size());
    memcpy(tree._version_roots[i].data(), _version_roots[i].data(),
      sizeof(HypernodeID) * _version_roots[i].size());
  }
  tree._out_degrees.resize(_out_degrees.size());
  for ( size_t i = 0; i < _out_degrees.size(); ++i ) {
    tree._out_degrees[i] = _out_degrees[i];
  }
  tree._incidence_array.resize(_incidence_array.size());
  memcpy(tree._incidence_array.data(), _incidence_array.data(),
    sizeof(HypernodeID) * _incidence_array.size());

  return tree;
}

// ! Resets internal data structures
void ContractionTree::reset() {
  tbb::parallel_invoke([&] {
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID hn) {
      _tree[hn].reset(hn);
      _out_degrees[hn].store(0);
    });
    _out_degrees[_num_hypernodes].store(0);
  }, [&] {
    parallel::parallel_free(_version_roots);
    _roots.clear();
  });
  _finalized = false;
}

// ! Free internal data in parallel
void ContractionTree::freeInternalData() {
  if ( _num_hypernodes > 0 ) {
    parallel::parallel_free(_tree, _roots, _out_degrees, _incidence_array);
  }
  _num_hypernodes = 0;
  _finalized = false;
}

void ContractionTree::memoryConsumption(utils::MemoryTreeNode* parent) const {
  ASSERT(parent);

  parent->addChild("Tree Nodes", sizeof(Node) * _tree.size());
  parent->addChild("Roots", sizeof(HypernodeID) * _roots.size());
  parent->addChild("Out-Degrees", sizeof(HypernodeID) * _out_degrees.size());
  parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
}

}  // namespace ds
}  // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2024 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <boost/core/span.hpp>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/static_graph.h"
#include "kahypar-resources/datastructure/fast_reset_flag_array.h"

using namespace mt_kahypar;
using FastResetArray = kahypar::ds::FastResetFlagArray<>;

using StaticGraph = ds::StaticGraph;

struct NeighborhoodResult {
  std::array<HypernodeID, 2> roots;
  const std::vector<HypernodeID>& n1_list;
  const FastResetArray& n1_set;
  const std::vector<HypernodeID>& n2_list;
  const FastResetArray& n2_set;
  bool includes_two_hop;

  bool isRoot(HypernodeID node) const {
    return node == roots[0] || node == roots[1];
  }

  bool isInN1(HypernodeID node) const {
    return isRoot(node) || n1_set[node];
  }

  bool isInN1Exactly(HypernodeID node) const {
    return n1_set[node];
  }

  bool isInN2(HypernodeID node) const {
    return isInN1(node) || n2_set[node];
  }

  bool isInN2Exactly(HypernodeID node) const {
    return n2_set[node];
  }
};


class NeighborhoodComputation {
 public:
  NeighborhoodComputation(HypernodeID num_nodes) {
    n1_set.setSize(num_nodes);
    n2_set.setSize(num_nodes);
  }

  template<size_t N>
  NeighborhoodResult computeNeighborhood(const StaticGraph& graph, std::array<HypernodeID, N> roots, bool include_two_hop) {
    return computeNeighborhood(graph, roots, include_two_hop, [](HypernodeID){ return true; });
  }

  template<size_t N, typename F>
  NeighborhoodResult computeNeighborhood(const StaticGraph& graph, std::array<HypernodeID, N> roots, bool include_two_hop, F filter) {
    static_assert(N > 0 && N <= 2);
    ALWAYS_ASSERT(n1_list.empty());
    NeighborhoodResult result {{roots[0], roots[0]}, n1_list, n1_set, n2_list, n2_set, include_two_hop};
    if constexpr (N == 2) {
      result.roots = roots;
    }

    for (HypernodeID root: roots) {
      for (HyperedgeID edge: graph.incidentEdges(root)) {
        HypernodeID neighbor = graph.edgeTarget(edge);
        if (!result.isInN1(neighbor) && filter(neighbor)) {
          n1_list.push_back(neighbor);
          n1_set.set(neighbor);
        }
      }
    }
    if (include_two_hop) {
      for (HypernodeID node: n1_list) {
        for (HyperedgeID edge: graph.incidentEdges(node)) {
          HypernodeID neighbor = graph.edgeTarget(edge);
          if (!result.isInN2(neighbor) && filter(neighbor)) {
            n2_list.push_back(neighbor);
            n2_set.set(neighbor);
          }
        }
      }
    }
    return result;
  }

  void reset() {
    n1_list.clear();
    n1_set.reset();
    n2_list.clear();
    n2_set.reset();
  }

 private:
  std::vector<HypernodeID> n1_list;
  FastResetArray n1_set;
  std::vector<HypernodeID> n2_list;
  FastResetArray n2_set;
};

class CliqueComputation {
 public:
  CliqueComputation(HypernodeID num_nodes) {
    local_neighbors.resize(num_nodes);
    current_set.setSize(num_nodes);
    tmp_set.setSize(num_nodes);
    forbidden.setSize(num_nodes);
    child = nullptr;
  }

  uint64_t computeMaxCliqueSize(const StaticGraph& graph, const std::vector<HypernodeID>& nodes) {
    current_set.reset();
    forbidden.reset();
    list.resize(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
      list[i] = nodes[i];
    }
    for (HypernodeID node: list) {
      current_set.set(node);
    }
    for (HypernodeID node: list) {
      local_neighbors[node].clear();
      for (HyperedgeID edge: graph.incidentEdges(node)) {
        HypernodeID neighbor = graph.edgeTarget(edge);
        if (current_set[neighbor]) {
          local_neighbors[node].push_back(neighbor);
        }
      }
    }
    return maxCliqueSizeRecursive<2>(nodes.size(), nodes.size(), 0);
  }

 private:
  template<size_t BRANCHING_FACTOR>
  uint64_t maxCliqueSizeRecursive(size_t original_nodes, size_t num_rem, size_t depth) {
    if (num_rem == 0) {
      return 0;
    } else if (num_rem == 1) {
      return 1;
    }
    else if ((depth % 2) == 0 && num_rem < original_nodes / 5 && num_rem >= 50) {
      if (child == nullptr) {
        child = std::make_unique<CliqueComputation>(local_neighbors.size());
      }
      auto r_nodes = remNodes(num_rem);
      ALWAYS_ASSERT(r_nodes.size() == num_rem);
      child->current_set.reset();
      child->forbidden.reset();
      child->list.resize(num_rem);
      for (size_t i = 0; i < r_nodes.size(); ++i) {
        child->list[i] = r_nodes[i];
      }
      for (HypernodeID node: child->list) {
        child->current_set.set(node);
      }
      for (HypernodeID node: child->list) {
        child->local_neighbors[node].clear();
        for (HypernodeID neighbor: local_neighbors[node]) {
          if (child->current_set[neighbor]) {
            child->local_neighbors[node].push_back(neighbor);
          }
        }
      }
      return child->maxCliqueSizeRecursive<BRANCHING_FACTOR>(num_rem, num_rem, depth);
    }

    std::array<std::pair<HypernodeID, HypernodeID>, BRANCHING_FACTOR> candidates;
    for (size_t i = 0; i < BRANCHING_FACTOR; ++i) {
      candidates[i] = {kInvalidHypernode, 0};
    }
    // compute prioritized candidates
    for (HypernodeID node: remNodes(num_rem)) {
      if (forbidden[node]) {
        // do some symmetry breaking
        continue;
      }

      HypernodeID count = 1;
      for (HypernodeID neighbor: local_neighbors[node]) {
        if (current_set[neighbor] && !forbidden[neighbor]) {
          count += 1;
        }
      }
      for (size_t i = 0; i < BRANCHING_FACTOR; ++i) {
        if (count > candidates[i].second) {
          auto tmp = std::make_pair(node, count);
          for (size_t j = i; j < BRANCHING_FACTOR; ++j) {
            std::swap(tmp, candidates[j]);
          }
          break;
        }
      }
    }

    uint64_t best = 0;
    for (auto [node, count]: candidates) {
      if (node == kInvalidHypernode) {
        break;
      }
      if (count > best) {
        // prepare data for recursion
        tmp_set.reset();
        for (HypernodeID neighbor: local_neighbors[node]) {
          tmp_set.set(neighbor);
        }
        auto node_list = remNodes(num_rem);
        auto it_after = std::partition(node_list.begin(), node_list.end(),
                          [&](HypernodeID curr_node) { return tmp_set[curr_node] && !forbidden[curr_node]; });
        size_t new_remaining = it_after - node_list.begin();
        if (new_remaining + 1 > best) {
          for (HypernodeID excluded: boost::span<HypernodeID>{it_after, node_list.end()}) {
            current_set.set(excluded, false);
          }
          // recursive call
          uint64_t new_val;
          if (depth < 3) {
            new_val = maxCliqueSizeRecursive<2>(original_nodes, new_remaining, depth + 1) + 1;
          } else {
            new_val = maxCliqueSizeRecursive<1>(original_nodes, new_remaining, depth + 1) + 1;
          }
          best = std::max(best, new_val);
          // restore data
          for (HypernodeID excluded: boost::span<HypernodeID>{it_after, node_list.end()}) {
            current_set.set(excluded, true);
          }
        }
      }
      // symmetry breaking for next candidates
      forbidden.set(node, true);
    }
    for (auto [node, count]: candidates) {
      if (node == kInvalidHypernode) {
        break;
      }
      forbidden.set(node, false);
    }
    return best;
  }

  boost::span<HypernodeID> remNodes(size_t num_rem) {
    ALWAYS_ASSERT(num_rem <= list.size());
    return boost::span<HypernodeID>{list.data(), num_rem};
  }

  std::vector<HypernodeID> list;
  std::vector<std::vector<HypernodeID>> local_neighbors;
  FastResetArray current_set;
  FastResetArray tmp_set;
  FastResetArray forbidden;
  std::unique_ptr<CliqueComputation> child;
};

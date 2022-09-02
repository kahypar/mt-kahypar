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

#include "snodes_coarsening_pass.h"


#include <algorithm>

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace star_partitioning {
SNodesCoarseningStage previous(SNodesCoarseningStage stage) {
  if (static_cast<uint8_t>(stage) > 0 && stage != SNodesCoarseningStage::ON_LARGE_GRAPH) {
    return static_cast<SNodesCoarseningStage>(static_cast<uint8_t>(stage) - 1);
  }
  return stage;
}

struct SNodesCoarseningPass::EqualityHash {
  static constexpr bool requires_check = false;
  using HashResult = HashValue;

  explicit EqualityHash() {
    seed = utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu());
  }

  HashResult calculateHash(vec<HypernodeID>& nodes) const {
    HashFunc hash(seed);
    std::sort(nodes.begin(), nodes.end());
    HashValue hash_value = kEdgeHashSeed;
    for ( const HypernodeID& node: nodes ) {
      HashValue next_hash = hash(node);
      hash_value *= 37;
      hash_value += 5;
      hash_value ^= next_hash;
    }
    return hash_value;
  }

  size_t combineHash(const HashResult& hash) const {
    return hash;
  }

 private:
  uint32_t seed;
};

// policies for different stages
bool appliesDegreeTwo(const SNodesCoarseningStage& stage) {
  return (stage != SNodesCoarseningStage::ANYTHING)
      && (stage != SNodesCoarseningStage::ON_LARGE_GRAPH)
      && static_cast<uint8_t>(stage) >= 1;
}

uint32_t treeDistanceForStage(const SNodesCoarseningStage& stage, bool use_spanning_tree) {
  if (use_spanning_tree && static_cast<uint8_t>(stage) <= 1) {
    return 3;
  } else if (use_spanning_tree && static_cast<uint8_t>(stage) <= 2) {
    return 6;
  }
  return 0;
}

bool appliesHighDegree(const SNodesCoarseningStage& stage) {
  return (stage != SNodesCoarseningStage::ANYTHING)
      && (stage != SNodesCoarseningStage::ON_LARGE_GRAPH)
      && static_cast<uint8_t>(stage) >= 4;
}

bool appliesTwins(const SNodesCoarseningStage& stage) {
  return (stage != SNodesCoarseningStage::ANYTHING);
}

bool appliesSimilarity(const SNodesCoarseningStage& stage) {
  return (stage != SNodesCoarseningStage::ANYTHING) && static_cast<uint8_t>(stage) >= 2;
}

bool relaxedSimilartiy(const SNodesCoarseningStage& stage) {
  return (stage != SNodesCoarseningStage::ANYTHING)
      && (stage != SNodesCoarseningStage::ON_LARGE_GRAPH)
      && static_cast<uint8_t>(stage) >= 3;
}

bool mixesDegreeOneAndTwo(const SNodesCoarseningStage& stage) {
  return (stage != SNodesCoarseningStage::ANYTHING)
      && (stage != SNodesCoarseningStage::ON_LARGE_GRAPH)
      && static_cast<uint8_t>(stage) >= 3;
}

double densityDiffForStage(const SNodesCoarseningStage& stage) {
  static const double diffs[7] = {1.4, 1.8, 2.1, 2.6, 4.1, 1000, 1.8};
  return diffs[static_cast<uint8_t>(stage)];
}

SNodesCoarseningPass::SNodesCoarseningPass(const Hypergraph& hg, const Context& context,
                                           const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage) :
  _hg(hg),
  _context(context),
  _s_nodes(_hg.separatedNodes().coarsest()),
  _part_ids(nullptr),
  _node_info_begin(),
  _node_info(),
  _current_num_nodes(_s_nodes.numVisibleNodes()),
  _target_num_nodes(target_num_nodes),
  _stage(stage) {
    ASSERT(_hg.initialNumNodes() == _s_nodes.numGraphNodes());
    tbb::parallel_invoke([&] {
      _node_info_begin.assign(_hg.initialNumNodes() + 1, 0);
    }, [&] {
      _node_info.assign(_current_num_nodes, NodeInfo());
    });
  }

SNodesCoarseningPass::SNodesCoarseningPass(const SeparatedNodes& s_nodes, const Hypergraph& hg, const Context& context,
                                           const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage) :
  _hg(hg),
  _context(context),
  _s_nodes(s_nodes),
  _part_ids(nullptr),
  _node_info_begin(),
  _node_info(),
  _current_num_nodes(_s_nodes.numVisibleNodes()),
  _target_num_nodes(target_num_nodes),
  _stage(stage) {
    ASSERT(_hg.initialNumNodes() == _s_nodes.numGraphNodes());
    tbb::parallel_invoke([&] {
      _node_info_begin.assign(_hg.initialNumNodes() + 1, 0);
    }, [&] {
      _node_info.assign(_current_num_nodes, NodeInfo());
    });
  }

HypernodeID SNodesCoarseningPass::run(vec<HypernodeID>& communities, bool use_spanning_tree) {
  tbb::parallel_invoke([&] {
    communities.clear();
    communities.assign(_current_num_nodes, kInvalidHypernode);
  }, [&] {
    setupNodeInfo(use_spanning_tree && _hg.initialNumNodes() < SpanningTree::max_num_nodes);
  });

  HypernodeID num_matches = runCurrentStage(communities, use_spanning_tree, true);
  while (_current_num_nodes - num_matches > _target_num_nodes
         && static_cast<uint8_t>(_stage) < static_cast<uint8_t>(SNodesCoarseningStage::ANYTHING)) {
    _stage = static_cast<SNodesCoarseningStage>(static_cast<uint8_t>(_stage) + 1);
    num_matches += runCurrentStage(communities, use_spanning_tree);
  }

  // assign any unmatched node to itself
  tbb::parallel_for(0UL, communities.size(), [&](const size_t i) {
    if (communities[i] == kInvalidHypernode) {
      communities[i] = i;
    }
  });

  ASSERT([&] {
    HypernodeID count = 0;
    for (HypernodeID i = 0; i < communities.size(); ++i) {
      if (communities[i] != i) {
        count++;
      }
    }
    return num_matches == count;
  }());
  return num_matches;
}

void SNodesCoarseningPass::setupNodeInfo(bool use_spanning_tree) {
  ASSERT(_s_nodes.outwardEdgesInitialized());

  Array<NodeInfo> tmp_node_info;
  tmp_node_info.assign(_s_nodes.numVisibleNodes(), NodeInfo());
  Array<parallel::IntegralAtomicWrapper<HypernodeID>> num_assigned_nodes;
  // sentinel for counting up later
  num_assigned_nodes.assign(_hg.initialNumNodes() + 1, parallel::IntegralAtomicWrapper<HypernodeID>(0));
  tbb::parallel_for(ID(0), _s_nodes.numVisibleNodes(), [&](const HypernodeID node) {
    HyperedgeWeight incident_weight_sum = 0;
    HyperedgeWeight max_weight = 0;
    HypernodeID max_target = kInvalidHypernode;
    for (const SeparatedNodes::Edge& e: _s_nodes.inwardEdges(node)) {
      incident_weight_sum += e.weight;
      if (e.weight >= max_weight) {
        max_weight = e.weight;
        max_target = e.target;
      }
    }
    NodeInfo& info = tmp_node_info[node];
    info.node = node;
    info.degree = _s_nodes.inwardDegree(node);
    info.density = _s_nodes.nodeWeight(node) == 0 ? std::numeric_limits<double>::infinity()
                        : static_cast<double>(incident_weight_sum) / static_cast<double>(_s_nodes.nodeWeight(node));
    info.assigned_graph_node = max_target;
    if (max_target != kInvalidHypernode) {
      num_assigned_nodes[max_target + 1].fetch_add(1);
    }
  });

  parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<HypernodeID>, Array>
          num_assigned_prefix_sum(num_assigned_nodes);

  tbb::parallel_invoke([&] {
    tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), _hg.initialNumNodes() + 1), num_assigned_prefix_sum);
    tbb::parallel_for(ID(0), _hg.initialNumNodes() + 1, [&](const HypernodeID node) {
      _node_info_begin[node] = num_assigned_nodes[node].load();
    });
  }, [&] {
    // Note: If the graph is too large, we don't construct a spanning tree. In this case,
    // the paths are initialized to kInvalidPath, but we also allow arbitrary distance.
    // This means that degree two nodes will be coarsened without regards to the tree distance.
    if (use_spanning_tree) {
      SpanningTree tree = constructMaxSpanningTree(_hg, TreePath::max_path_length);
      tree.calculatePaths([&](const HypernodeID& graph_node, TreePath path) {
        for (const SeparatedNodes::Edge& e: _s_nodes.outwardEdges(graph_node)) {
          if (e.target < _s_nodes.numVisibleNodes()) {
            NodeInfo& info = tmp_node_info[e.target];
            if (info.degree == 2 && info.assigned_graph_node != graph_node) {
              info.tree_path = path;
            }
          }
        }
      });
    }
  });

  tbb::parallel_for(ID(0), _s_nodes.numVisibleNodes(), [&](const HypernodeID node) {
    const NodeInfo& tmp_info = tmp_node_info[node];
    if (tmp_info.assigned_graph_node != kInvalidHypernode) {
      const HypernodeID new_index = num_assigned_nodes[tmp_info.assigned_graph_node].fetch_add(1);
      _node_info[new_index] = tmp_info;
    } else {
      const HypernodeID new_index = num_assigned_nodes[_hg.initialNumNodes()].fetch_add(1);
      _node_info[new_index] = tmp_info;
    }
  });
  ASSERT(_node_info.size() == _s_nodes.numVisibleNodes());
}

HypernodeID SNodesCoarseningPass::runCurrentStage(vec<HypernodeID>& communities, bool use_spanning_tree, bool first) {
  LocalizedData data;
  Params params;
  params.max_node_weight = _context.coarsening.max_allowed_node_weight;
  params.degree_one_cluster_size = (_stage == SNodesCoarseningStage::D1_TWINS) ? MAX_CLUSTER_SIZE : 2;
  params.accepted_density_diff = densityDiffForStage(_stage);
  params.max_tree_distance = treeDistanceForStage(_stage, use_spanning_tree);

  if (first) {
    applyDegreeZeroCoarsening(params, communities, data);
  }
  if (_stage != SNodesCoarseningStage::ANYTHING) {
    if (appliesTwins(_stage)) {
      applyHashingRound<EqualityHash>(params, communities, data, 2);
    }
    if (appliesSimilarity(_stage)) {
      for (size_t i = 0; i < NUM_SIMILARITY_ROUNDS; ++i) {
        if (relaxedSimilartiy(_stage)) {
          applyHashingRound<SimilarityHash<2>>(params, communities, data, 3, SIMILARITY_RATIO_RELAXED);
        } else {
          applyHashingRound<SimilarityHash<3>>(params, communities, data, 3, SIMILARITY_RATIO);
        }
      }
    }

    _hg.doParallelForAllNodes([&](const HypernodeID& node) {
      applyCoarseningForNode(params, communities, data, node);
    });
  } else {
    _hg.doParallelForAllNodes([&](const HypernodeID& node) {
      HyperedgeID& counter = data.match_counter.local();
      vec<HypernodeID>& any_nodes = data.high_degree_nodes.local();
      any_nodes.clear();
      for (HypernodeID index = _node_info_begin[node]; index < _node_info_begin[node + 1]; ++index) {
        ASSERT(info(index).assigned_graph_node == node && info(index).density > 0);
        if (communities[info(index).node] == kInvalidHypernode) {
          any_nodes.push_back(index);
        }
      }

      sortByDensity(any_nodes);
      matchPairs(params, communities, any_nodes, counter);
    });
  }
  return data.match_counter.combine(std::plus<>());
}

void SNodesCoarseningPass::applyDegreeZeroCoarsening(const Params& params, vec<HypernodeID>& communities, LocalizedData& data) {
  const HypernodeID first_d0 = _node_info_begin.back();
  const HypernodeID last_d0 = _node_info.size();
  if (_part_ids != nullptr) {
    std::sort(_node_info.begin() + first_d0, _node_info.end(), [&] (const NodeInfo& l, const NodeInfo& r) {
      return _part_ids[l.node] < _part_ids[r.node];
    });
  }
  const HypernodeID num_clusters = (last_d0 - first_d0 + DEGREE_ZERO_CLUSTER_SIZE - 1) / DEGREE_ZERO_CLUSTER_SIZE;
  if (last_d0 - first_d0 > 4) { // TODO: magic number
    tbb::parallel_for(ID(0), num_clusters, [&](const HypernodeID cluster) {
      const HypernodeID first_node_of_block = first_d0 + cluster * DEGREE_ZERO_CLUSTER_SIZE;
      HypernodeID offset = 0;
      while (offset < DEGREE_ZERO_CLUSTER_SIZE && first_node_of_block + offset < last_d0) {
        const HypernodeID start = first_node_of_block + offset;
        HypernodeID& counter = data.match_counter.local();
        HypernodeWeight weight = 0;
        for (HypernodeID i = 0; offset + i < DEGREE_ZERO_CLUSTER_SIZE && start + i < last_d0; ++i) {
          const HypernodeID node = info(start + i).node;
          if (weight + _s_nodes.nodeWeight(node) <= params.max_node_weight
              && (_part_ids == nullptr || _part_ids[start] == _part_ids[node])) {
            ++offset;
            ASSERT(communities[node] == kInvalidHypernode);
            communities[node] = info(start).node;
            weight += _s_nodes.nodeWeight(node);
            if (i > 0) {
              ++counter;
            }
          } else if (i > 0) {
            break;
          } else {
            ++offset;
            break;
          }
        }
      }
    });
  }
}

void SNodesCoarseningPass::applyCoarseningForNode(const Params& params, vec<HypernodeID>& communities,
                                                  LocalizedData& data, const HypernodeID& node) {
  HyperedgeID& counter = data.match_counter.local();
  vec<HypernodeID>& degree_one = data.degree_one_nodes.local();
  degree_one.clear();
  vec<HypernodeID>& degree_two = data.degree_two_nodes.local();
  degree_two.clear();
  vec<HypernodeID>& high_degree = data.high_degree_nodes.local();
  high_degree.clear();

  for (HypernodeID index = _node_info_begin[node]; index < _node_info_begin[node + 1]; ++index) {
    ASSERT(info(index).assigned_graph_node == node && info(index).density > 0);
    if (communities[info(index).node] == kInvalidHypernode) {
      if (info(index).degree == 1) {
        degree_one.push_back(index);
      } else if (info(index).degree == 2) {
        degree_two.push_back(index);
      } else {
        high_degree.push_back(index);
      }
    }
  }

  // degree one nodes
  sortByDensity(degree_one);
  size_t current_index = 0;
  while (current_index < degree_one.size()) {
    const double starting_density = info(degree_one[current_index]).density;
    const HypernodeID starting_node = info(degree_one[current_index]).node;
    HypernodeID cluster_size = 1;
    HypernodeWeight weight = 0;
    ++current_index;
    while (cluster_size < params.degree_one_cluster_size
           && current_index < degree_one.size()
           && info(degree_one[current_index]).density / starting_density <= params.accepted_density_diff
           && weight + _s_nodes.nodeWeight(info(degree_one[current_index]).node) <= params.max_node_weight
           && (_part_ids == nullptr || _part_ids[starting_node] == _part_ids[info(degree_one[current_index]).node])) {
      communities[info(degree_one[current_index]).node] = starting_node;
      communities[starting_node] = starting_node;
      weight += _s_nodes.nodeWeight(info(degree_one[current_index]).node);
      ++current_index;
      ++cluster_size;
    }
    counter += (cluster_size - 1);
  }

  if (appliesDegreeTwo(_stage)) {
    // degree two nodes
    std::sort(degree_two.begin(), degree_two.end(), [&](const HypernodeID& left, const HypernodeID& right) {
      if (_part_ids == nullptr || _part_ids[info(left).node] == _part_ids[info(right).node]) {
        if (info(left).tree_path == info(right).tree_path) {
          return info(left).density < info(right).density;
        }
        return info(left).tree_path < info(right).tree_path;
      }
      return _part_ids[info(left).node] < _part_ids[info(right).node];
    });
    for (size_t i = 0; i + 1 < degree_two.size(); ++i) {
      for (size_t j = 1; j <= D2_SEARCH_RANGE && i + j < degree_two.size(); ++j) {
        const NodeInfo& curr = info(degree_two[i]);
        const NodeInfo& next = info(degree_two[i + j]);
        if (std::max(curr.density / next.density, next.density / curr.density) <= params.accepted_density_diff
            && (params.max_tree_distance == 0 || curr.tree_path.distance(next.tree_path) <= params.max_tree_distance)
            && validPair(params, curr.node, next.node)) {
          ASSERT(communities[curr.node] == kInvalidHypernode && communities[next.node] == kInvalidHypernode);
          communities[curr.node] = curr.node;
          communities[next.node] = curr.node;
          std::swap(degree_two[i + j], degree_two[i + 1]);
          ++counter;
          ++i;
          break;
        }
      }
    }
  }

  if (appliesHighDegree(_stage)) {
    sortByDensity(high_degree);
    matchPairs(params, communities, high_degree, counter);
  }
}

void SNodesCoarseningPass::sortByDensity(vec<HypernodeID>& nodes) {
  std::sort(nodes.begin(), nodes.end(), [&](const HypernodeID& left, const HypernodeID& right) {
    if (_part_ids == nullptr || _part_ids[info(left).node] == _part_ids[info(right).node]) {
      return info(left).density < info(right).density;
    }
    return _part_ids[info(left).node] < _part_ids[info(right).node];
  });
}

void SNodesCoarseningPass::matchPairs(const Params& params, vec<HypernodeID>& communities,
                                      vec<HypernodeID>& nodes, HypernodeID& counter) {
  for (size_t i = 0; i + 1 < nodes.size(); ++i) {
    const NodeInfo& curr = info(nodes[i]);
    const NodeInfo& next = info(nodes[i + 1]);
    if (std::max(curr.density / next.density, next.density / curr.density) <= params.accepted_density_diff
        && validPair(params, curr.node, next.node)) {
      ASSERT(communities[curr.node] == kInvalidHypernode && communities[next.node] == kInvalidHypernode);
      communities[curr.node] = curr.node;
      communities[next.node] = curr.node;
      ++counter;
      ++i;
    }
  }
}

bool SNodesCoarseningPass::validPair(const Params& params, const HypernodeID& l, const HypernodeID& r) {
  if (_s_nodes.nodeWeight(l) + _s_nodes.nodeWeight(r) > params.max_node_weight) {
    return false;
  }
  return _part_ids == nullptr || _part_ids[l] == _part_ids[r];
}

std::pair<HyperedgeWeight, HyperedgeWeight> SNodesCoarseningPass::intersection_and_union(
            const HypernodeID& s_node_left, const HypernodeID& s_node_right,
            vec<std::pair<HypernodeID, HyperedgeWeight>>& lhs,
            vec<std::pair<HypernodeID, HyperedgeWeight>>& rhs) {
  lhs.clear();
  for (const SeparatedNodes::Edge& e: _s_nodes.inwardEdges(s_node_left)) {
    lhs.emplace_back(e.target, e.weight);
  }
  rhs.clear();
  for (const SeparatedNodes::Edge& e: _s_nodes.inwardEdges(s_node_right)) {
    rhs.emplace_back(e.target, e.weight);
  }
  return intersection_and_union(lhs, rhs);
}

std::pair<HyperedgeWeight, HyperedgeWeight>
SNodesCoarseningPass::intersection_and_union(vec<std::pair<HypernodeID, HyperedgeWeight>>& lhs,
                                             vec<std::pair<HypernodeID, HyperedgeWeight>>& rhs) {
  std::sort(lhs.begin(), lhs.end(), [](const auto& l, const auto& r) {
    return l.first < r.first;
  });
  lhs.emplace_back(std::numeric_limits<HypernodeID>::max(), 0);
  std::sort(rhs.begin(), rhs.end(), [](const auto& l, const auto& r) {
    return l.first < r.first;
  });
  rhs.emplace_back(std::numeric_limits<HypernodeID>::max(), 0);
  HyperedgeWeight intersection_weight = 0;
  HyperedgeWeight union_weight = 0;
  size_t i = 0;
  size_t j = 0;
  while ( i < lhs.size() && j < rhs.size() ) {
    if ( lhs[i].first == rhs[j].first ) {
      union_weight += std::max(lhs[i].second, rhs[j].second);
      intersection_weight += std::min(lhs[i].second, rhs[j].second);
      ++i;
      ++j;
    } else if ( lhs[i].first < rhs[j].first ) {
      union_weight += lhs[i].second;
      ++i;
    } else {
      union_weight += rhs[j].second;
      ++j;
    }
  }
  return {intersection_weight, union_weight};
}

} // namepace star_partitioning
} // namespace mt_kahypar

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
using HashFunc = kahypar::math::MurmurHash<HypernodeID>;
using HashValue = typename HashFunc::HashValue;

struct SNodesCoarseningPass::Footprint {
  bool operator==(const Footprint& other) const {
    for ( size_t i = 0; i < NUM_HASHES; ++i ) {
      if ( footprint[i] != other.footprint[i] ) {
        return false;
      }
    }
    return true;
  }

  bool operator<(const Footprint& other) const {
    for ( size_t i = 0; i < NUM_HASHES; ++i ) {
      if ( footprint[i] < other.footprint[i] ) {
        return true;
      } else if ( footprint[i] > other.footprint[i] ) {
        return false;
      }
    }
    return false;
  }

  HashValue footprint[NUM_HASHES];
};

// a similarity hashing function (see also hypergraph_sparsifier.h)
struct SNodesCoarseningPass::SimilarityHash {
  static constexpr bool requires_check = true;
  using HashResult = Footprint;

  explicit SimilarityHash() {
    for (size_t i = 0; i < NUM_HASHES; ++i) {
      int seed = utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu());
      hash_functions[i] = HashFunc(seed);
    }
  }

  HashResult calculateHash(const vec<HypernodeID>& nodes) const {
    Footprint footprint;
    for ( size_t i = 0; i < NUM_HASHES; ++i ) {
      footprint.footprint[i] = minHash(hash_functions[i], nodes);
    }
    return footprint;
  }

  size_t combineHash(const HashResult& footprint) const {
    HashValue hash_value = kEdgeHashSeed;
    for ( const HashValue& value : footprint.footprint ) {
      hash_value ^= value;
    }
    return hash_value;
  }

 private:
  HashValue minHash(const HashFunc& hash_function,
                    const vec<HypernodeID>& nodes ) const {
    HashValue hash_value = std::numeric_limits<HashValue>::max();
    for ( const HypernodeID& node : nodes ) {
      hash_value = std::min(hash_value, hash_function(node));
    }
    return hash_value;
  }

  HashFunc hash_functions[NUM_HASHES];
};

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

bool allowsDegreeTwo(const SNodesCoarseningStage& stage) {
  return static_cast<uint8_t>(stage) >= 3;
}

bool allowsHighDegree(const SNodesCoarseningStage& stage) {
  return static_cast<uint8_t>(stage) >= 4;
}

bool appliesTwins(const SNodesCoarseningStage& stage) {
  // Note: We don't apply twins at later stages because most probably there
  // aren't any twins left anymore.
  return stage == SNodesCoarseningStage::DEGREE_ONE_AND_TWINS
      || stage == SNodesCoarseningStage::DEGREE_ONE_AND_DEGREE_TWO_AND_TWINS;
}

SNodesCoarseningPass::SNodesCoarseningPass(const Hypergraph& hg, const Context& context,
                                           const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage) :
  _hg(hg),
  _context(context),
  _s_nodes(_hg.separatedNodes().onliest()),
  _node_info_begin(),
  _node_info(),
  _current_num_nodes(_s_nodes.numNodes()),
  _target_num_nodes(target_num_nodes),
  _stage(stage) {
    tbb::parallel_invoke([&] {
      _node_info_begin.assign(_hg.initialNumNodes() + 1, 0);
    }, [&] {
      _node_info.assign(_current_num_nodes, FullNodeInfo());
    });
  }

void SNodesCoarseningPass::run(vec<HypernodeID>& communities) {
  tbb::parallel_invoke([&] {
    communities.clear();
    communities.assign(_current_num_nodes, kInvalidHypernode);
  }, [&] {
    setupNodeInfo();
  });

  HypernodeID num_matches = runCurrentStage(communities, true);
  while (_current_num_nodes - num_matches > _target_num_nodes
         && _stage != SNodesCoarseningStage::ANYTHING) {
    _stage = static_cast<SNodesCoarseningStage>(static_cast<uint8_t>(_stage) + 1);
    num_matches += runCurrentStage(communities);
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
}

void SNodesCoarseningPass::setupNodeInfo() {
  Array<NodeInfo> tmp_node_info;
  tmp_node_info.assign(_s_nodes.numNodes(), NodeInfo());
  Array<parallel::IntegralAtomicWrapper<HypernodeID>> num_assigned_nodes;
  // sentinel for counting up later
  num_assigned_nodes.assign(_hg.initialNumNodes() + 1, parallel::IntegralAtomicWrapper<HypernodeID>(0));
  tbb::parallel_for(ID(0), _s_nodes.numNodes(), [&](const HypernodeID node) {
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
    info.density = _s_nodes.nodeWeight(node) == 0 ? std::numeric_limits<double>::infinity()
                        : static_cast<double>(incident_weight_sum) / static_cast<double>(_s_nodes.nodeWeight(node));
    info.assigned_graph_node = max_target;
    if (max_target != kInvalidHypernode) {
      num_assigned_nodes[max_target + 1].fetch_add(1);
    }
  });

  parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<HypernodeID>, Array>
          num_assigned_prefix_sum(num_assigned_nodes);
  tbb::parallel_scan(tbb::blocked_range<size_t>(ID(0), _hg.initialNumNodes() + 1), num_assigned_prefix_sum);

  tbb::parallel_for(ID(0), _hg.initialNumNodes() + 1, [&](const HypernodeID node) {
    _node_info_begin[node] = num_assigned_nodes[node].load();
  });

  tbb::parallel_for(ID(0), _s_nodes.numNodes(), [&](const HypernodeID node) {
    const NodeInfo& tmp_info = tmp_node_info[node];
    if (tmp_info.assigned_graph_node != kInvalidHypernode) {
      const HypernodeID new_index = num_assigned_nodes[tmp_info.assigned_graph_node].fetch_add(1);
      _node_info[new_index] = FullNodeInfo(node, _s_nodes.inwardDegree(node), tmp_info);
    } else {
      const HypernodeID new_index = num_assigned_nodes[_hg.initialNumNodes()].fetch_add(1);
      _node_info[new_index] = FullNodeInfo(node, 0, tmp_info);
    }
  });
  ASSERT(_node_info.size() == _s_nodes.numNodes());
}

HypernodeID SNodesCoarseningPass::runCurrentStage(vec<HypernodeID>& communities, bool first) {
  LocalizedData data;
  Params params;
  params.max_node_weight = _context.coarsening.max_allowed_node_weight;
  params.degree_one_cluster_size = _stage == SNodesCoarseningStage::PREFERABLE_DEGREE_ONE ? MAX_CLUSTER_SIZE : 2;
  if (_stage == SNodesCoarseningStage::PREFERABLE_DEGREE_ONE) {
    params.accepted_density_diff = PREFERRED_DENSITY_DIFF;
  } else if (_stage == SNodesCoarseningStage::ANY_DEGREE_RELAXED) {
    params.accepted_density_diff = RELAXED_DENSITY_DIFF;
  } else if (_stage == SNodesCoarseningStage::ANYTHING) {
    params.accepted_density_diff = std::numeric_limits<double>::infinity();
  } else {
    params.accepted_density_diff = TOLERABLE_DENSITY_DIFF;
  }

  if (first) {
    applyDegreeZeroCoarsening(params, communities, data);
  }
  if (appliesTwins(_stage)) {
    applyHashingRound<EqualityHash>(params, communities, data, 2);
  }
  _hg.doParallelForAllNodes([&](const HypernodeID& node) {
    applyCoarseningForNode(params, communities, data, node);
  });
  return data.match_counter.combine(std::plus<>());
}

void SNodesCoarseningPass::applyDegreeZeroCoarsening(const Params& params, vec<HypernodeID>& communities, LocalizedData& data) {
  const HypernodeID first_d0 = _node_info_begin.back();
  const HypernodeID last_d0 = _node_info.size();
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
          if (weight + _s_nodes.nodeWeight(node) <= params.max_node_weight) {
            ++offset;
            ASSERT(communities[node] == kInvalidHypernode);
            communities[node] = info(start).node;
            weight += _s_nodes.nodeWeight(node);
            if (i > 0) {
              ++counter;
            }
          } else if (i > 0) {
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

  for (HypernodeID index = _node_info_begin[node]; index < _node_info_begin[node + 1]; ++index) {
    ASSERT(info(index).assigned_graph_node == node && info(index).density > 0);
    if (communities[info(index).node] == kInvalidHypernode) {
      if (info(index).degree == 1) {
        degree_one.push_back(index);
      } else if (info(index).degree == 2) {
        degree_two.push_back(index);
      } else {
        // TODO
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
           && weight + _s_nodes.nodeWeight(info(degree_one[current_index]).node) <= params.max_node_weight) {
      communities[info(degree_one[current_index]).node] = starting_node;
      communities[starting_node] = starting_node;
      weight += _s_nodes.nodeWeight(info(degree_one[current_index]).node);
      ++current_index;
      ++cluster_size;
    }
    counter += (cluster_size - 1);
  }
}

void SNodesCoarseningPass::sortByDensity(vec<HypernodeID>& nodes) {
  std::sort(nodes.begin(), nodes.end(), [&](const HypernodeID& left, const HypernodeID& right) {
    return info(left).density < info(right).density;
  });
}


std::pair<HyperedgeID, HyperedgeID> SNodesCoarseningPass::intersection_and_union(
            const HypernodeID& s_node_left, const HypernodeID& s_node_right,
            vec<HypernodeID>& lhs, vec<HypernodeID>& rhs) {
  lhs.clear();
  for (const SeparatedNodes::Edge& e: _s_nodes.inwardEdges(s_node_left)) {
    lhs.push_back(e.target);
  }
  std::sort(lhs.begin(), lhs.end());
  rhs.clear();
  for (const SeparatedNodes::Edge& e: _s_nodes.inwardEdges(s_node_right)) {
    rhs.push_back(e.target);
  }
  std::sort(rhs.begin(), rhs.end());
  HyperedgeID intersection_size = 0;
  size_t i = 0;
  size_t j = 0;
  while ( i < lhs.size() && j < rhs.size() ) {
    if ( lhs[i] == rhs[j] ) {
      ++intersection_size;
      ++i;
      ++j;
    } else if ( lhs[i] < rhs[j] ) {
      ++i;
    } else {
      ++j;
    }
  }
  const HyperedgeID union_size = lhs.size() + rhs.size() - intersection_size;
  return {intersection_size, union_size};
}

} // namepace star_partitioning
} // namespace mt_kahypar

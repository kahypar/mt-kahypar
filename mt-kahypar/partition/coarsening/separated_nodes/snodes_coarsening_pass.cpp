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

#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

namespace mt_kahypar {
namespace star_partitioning {

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
  // setup data
  tbb::parallel_invoke([&] {
    communities.clear();
    communities.assign(_current_num_nodes, kInvalidHypernode);
  }, [&] {
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
  });

}

} // namepace star_partitioning
} // namespace mt_kahypar

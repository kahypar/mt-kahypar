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

#include "gmock/gmock.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/snodes_coarsening_pass.h"

using ::testing::Test;

namespace mt_kahypar {
using star_partitioning::SNodesCoarseningPass;
using star_partitioning::SNodesCoarseningStage;
using ds::SeparatedNodes;
using ds::SepNodesStack;

class ACoarseningPass : public Test {
 template<typename T>
 using vec = parallel::scalable_vector<T>;

 public:
  ACoarseningPass() : graph(), context(), stack(), communities() { }

  void initialize(HypernodeID num_nodes, HyperedgeID num_edges, vec<vec<HypernodeID>> edges,
                  vec<vec<std::pair<HypernodeID, HyperedgeWeight>>> sep_nodes) {
    graph = HypergraphFactory::construct(num_nodes, num_edges, edges);
    stack = SepNodesStack(num_nodes);
    graph.setSeparatedNodes(&stack);

    SeparatedNodes& s_nodes = stack.onliest();
    vec<std::tuple<HypernodeID, HyperedgeID, HypernodeWeight>> nodes;
    vec<SeparatedNodes::Edge> s_edges;
    for (const auto& edges_of_node: sep_nodes) {
      nodes.emplace_back(kInvalidHypernode, s_edges.size(), 1);
      for (const auto& [target, weight]: edges_of_node) {
        s_edges.emplace_back(target, weight);
      }
    }
    s_nodes.addNodes(nodes, s_edges);
  }

  SNodesCoarseningPass setupPass(const HypernodeID& target_num_nodes, const SNodesCoarseningStage& stage) {
    return SNodesCoarseningPass(graph, context, target_num_nodes, stage);
  }

  Hypergraph graph;
  Context context;
  SepNodesStack stack;
  vec<HypernodeID> communities;
};

TEST_F(ACoarseningPass, setupData) {
  initialize(3, 2, {{0, 1}, {1, 2}}, { {}, {}, {{0, 1}}, {{0, 2}},
                                       {{0, 1}, {1, 2}}, {{1, 2}, {2, 1}}, {{0, 1}, {1, 1}, {2, 2}} });
  SNodesCoarseningPass c_pass = setupPass(4, SNodesCoarseningStage::DEGREE_ZERO);
  c_pass.run(communities);

  const auto& info = c_pass.nodeInfo();
  ASSERT_EQ(0, info[0].assigned_graph_node);
  ASSERT_EQ(0, info[1].assigned_graph_node);
  ASSERT_EQ(1, info[2].assigned_graph_node);
  ASSERT_EQ(1, info[3].assigned_graph_node);
  ASSERT_EQ(2, info[4].assigned_graph_node);
  ASSERT_EQ(kInvalidHypernode, info[5].assigned_graph_node);
  ASSERT_EQ(kInvalidHypernode, info[6].assigned_graph_node);
  ASSERT_EQ(3, info[2].density);
  ASSERT_EQ(3, info[3].density);
  ASSERT_EQ(4, info[4].density);
  ASSERT_EQ(0, info[5].density);
  ASSERT_EQ(0, info[6].density);

  const auto& info_begin = c_pass.nodeInfoBegin();
  ASSERT_EQ(0, info_begin[0]);
  ASSERT_EQ(2, info_begin[1]);
  ASSERT_EQ(4, info_begin[2]);
  ASSERT_EQ(5, info_begin[3]);
}

}  // namespace mt_kahypar

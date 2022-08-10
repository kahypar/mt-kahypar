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
using star_partitioning::kInvalidPath;
using ds::SeparatedNodes;
using ds::SepNodesStack;

class ACoarseningPass : public Test {
 template<typename T>
 using vec = parallel::scalable_vector<T>;

 public:
  ACoarseningPass() : graph(), context(), stack(), communities() {
    context.coarsening.max_allowed_node_weight = 10;
  }

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
    s_nodes.initializeOutwardEdges();
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
  SNodesCoarseningPass c_pass = setupPass(4, SNodesCoarseningStage::D1_TWINS);
  c_pass.run(communities);
  ASSERT_EQ(7, communities.size());

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

  ASSERT_EQ(kInvalidPath, info[0].tree_path);
  ASSERT_EQ(kInvalidPath, info[1].tree_path);
  ASSERT_NE(kInvalidPath, info[2].tree_path);
  ASSERT_NE(kInvalidPath, info[3].tree_path);
  ASSERT_EQ(kInvalidPath, info[4].tree_path);
  ASSERT_EQ(kInvalidPath, info[5].tree_path);
  ASSERT_EQ(kInvalidPath, info[6].tree_path);

  const auto& info_begin = c_pass.nodeInfoBegin();
  ASSERT_EQ(0, info_begin[0]);
  ASSERT_EQ(2, info_begin[1]);
  ASSERT_EQ(4, info_begin[2]);
  ASSERT_EQ(5, info_begin[3]);
}

TEST_F(ACoarseningPass, removesDegreeZero) {
  initialize(1, 0, {}, { {}, {}, {}, {}, {}, {} });
  context.coarsening.max_allowed_node_weight = 3;
  SNodesCoarseningPass c_pass = setupPass(2, SNodesCoarseningStage::D1_TWINS);
  c_pass.run(communities);

  HypernodeID c_0 = communities[0];
  HyperedgeID c_1 = kInvalidHypernode;
  HypernodeID num_c_0 = 1;
  for (size_t i = 1; i < communities.size(); ++i) {
    if (communities[i] == c_0) {
      ++num_c_0;
    } else if (c_1 == kInvalidHypernode) {
      c_1 = communities[i];
    } else {
      ASSERT_EQ(c_1, communities[i]);
    }
  }
  ASSERT_EQ(num_c_0, 3);
}

TEST_F(ACoarseningPass, coarsensDegreeOneNodes) {
  initialize(1, 0, {}, { {{0, 1}}, {{0, 2}}, {{0, 3}}, {{0, 6}}, {{0, 7}}, {{0, 8}}, {{0, 8}}, {{0, 9}} });
  SNodesCoarseningPass c_pass = setupPass(4, SNodesCoarseningStage::D1_TWINS);
  c_pass.run(communities);

  ASSERT_EQ(0, communities[0]);
  ASSERT_EQ(1, communities[1]);
  ASSERT_EQ(1, communities[2]);
  ASSERT_EQ(3, communities[3]);
  ASSERT_EQ(3, communities[4]);
  ASSERT_EQ(3, communities[5]);
  ASSERT_EQ(3, communities[6]);
  ASSERT_EQ(7, communities[7]);

  SNodesCoarseningPass c_pass_2 = setupPass(5, SNodesCoarseningStage::D1_D2_TWINS_SIMILARITY);
  c_pass_2.run(communities);

  ASSERT_EQ(0, communities[0]);
  ASSERT_EQ(0, communities[1]);
  ASSERT_EQ(2, communities[2]);
  ASSERT_EQ(2, communities[3]);
  ASSERT_EQ(4, communities[4]);
}

TEST_F(ACoarseningPass, findsTwins) {
  initialize(4, 0, {}, { {{0, 1}, {1, 1}}, {{0, 1}, {1, 1}},
                         {{0, 1}, {1, 1}, {2, 1}}, {{2, 1}, {1, 1}, {0, 1}}, {{2, 2}, {1, 1}, {0, 3}}, {{2, 2}, {0, 2}, {1, 1}},
                         {{0, 1}, {3, 1}, {2, 1}}, {{2, 1}, {3, 1}, {0, 1}}, {{1, 1}, {2, 1}, {3, 1}},
                         {{0, 1}, {1, 1}, {2, 1}, {3, 1}}, {{0, 3}, {1, 3}, {2, 3}, {3, 3}} });
  SNodesCoarseningPass c_pass_2 = setupPass(7, SNodesCoarseningStage::D1_D2_TWINS);
  c_pass_2.run(communities);

  ASSERT_EQ(communities[0], communities[1]);
  ASSERT_EQ(communities[2], communities[3]);
  ASSERT_EQ(communities[4], communities[5]);
  ASSERT_EQ(communities[6], communities[7]);
  ASSERT_EQ(8, communities[8]);
  ASSERT_EQ(9, communities[9]);
  ASSERT_EQ(10, communities[10]);
}

TEST_F(ACoarseningPass, coarsensDegreeTwoNodes) {
  initialize(6, 5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}},
             { {{0, 2}, {1, 1}}, {{0, 2}, {2, 1}}, {{0, 4}, {3, 3}}, {{0, 9}, {4, 8}}, {{0, 4}, {5, 3}} });
  SNodesCoarseningPass c_pass = setupPass(3, SNodesCoarseningStage::D1_D2_TWINS);
  c_pass.run(communities);

  ASSERT_EQ(communities[0], communities[1]);
  ASSERT_EQ(communities[2], communities[4]);
  ASSERT_EQ(communities[3], communities[3]);
}

}  // namespace mt_kahypar

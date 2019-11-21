/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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


#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/partition/initial_partitioning/kahypar.h"

using ::testing::Test;

namespace mt_kahypar {

class AKahyparPartitioner : public ds::AHypergraph<2> {

 private:
  using Base = ds::AHypergraph<2>;

 public:
  using Base::TBBArena;
  using Base::TestStreamingHypergraph;
  using Base::TestHypergraph;

  AKahyparPartitioner() :
    Base(),
    hypergraph(construct_hypergraph( 7,
      { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} },
      { 0, 0, 0, 1, 1, 1, 1 },
      { 0, 0, 1, 1 },
      { 0, 0, 1, 1, 2, 3, 2 } )),
    node_mapping({0, 1, 2, 3, 4, 5, 6}),
    reverse_mapping({0, 1, 2, 3, 4, 5, 6}) { }

  static void TearDownTestSuite() {
    TBBArena::instance().terminate();
  }

  TestHypergraph hypergraph;
  std::vector<HypernodeID> node_mapping;
  std::vector<HypernodeID> reverse_mapping;
};

template<typename T>
bool isEqual(const parallel::scalable_vector<T>& v1,
             const parallel::scalable_vector<T>& v2)
{
  return (v1.size() == v2.size() &&
          std::equal(v1.begin(), v1.end(), v2.begin()));
}


template< class HyperGraph >
void assignPartitionIDs(HyperGraph& hypergraph) {
  using StreamingHyperGraph = typename HyperGraph::StreamingHypergraph;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    PartitionID part_id = StreamingHyperGraph::get_numa_node_of_vertex(hn);
    hypergraph.setNodePart(hn, part_id);
  }
  hypergraph.updateGlobalPartInfos();
  hypergraph.initializeNumCutHyperedges();
}

TEST_F(AKahyparPartitioner, ConvertsHypergraphWithCorrectNumberOfVertices) {
  KaHyParHypergraph kahypar_hg = convertToKaHyParHypergraph(hypergraph, node_mapping);
  ASSERT_EQ(7, kahypar_hg.num_vertices);
}

TEST_F(AKahyparPartitioner, ConvertsHypergraphWithCorrectNumberOfHyperedges) {
  KaHyParHypergraph kahypar_hg = convertToKaHyParHypergraph(hypergraph, node_mapping);
  ASSERT_EQ(4, kahypar_hg.num_hyperedges);
}

TEST_F(AKahyparPartitioner, ConvertsHypergraphWithCorrectHyperedgeIndices) {
  KaHyParHypergraph kahypar_hg = convertToKaHyParHypergraph(hypergraph, node_mapping);
  ASSERT_TRUE(isEqual(kahypar_hg.hyperedge_indices, {0, 2, 6, 9, 12}));
}

TEST_F(AKahyparPartitioner, ConvertsHypergraphWithCorrectHyperedges) {
  KaHyParHypergraph kahypar_hg = convertToKaHyParHypergraph(hypergraph, node_mapping);
  ASSERT_TRUE(isEqual(kahypar_hg.hyperedges, {0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6}));
}

TEST_F(AKahyparPartitioner, ExtractsBlock0WithCorrectNumberOfVertices) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 0, node_mapping);
  ASSERT_EQ(3, kahypar_hg.num_vertices);
}

TEST_F(AKahyparPartitioner, ExtractsBlock1WithCorrectNumberOfVertices) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 1, node_mapping);
  ASSERT_EQ(4, kahypar_hg.num_vertices);
}

TEST_F(AKahyparPartitioner, ExtractsBlock0WithCorrectNumberOfHyperedges) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 0, node_mapping);
  ASSERT_EQ(3, kahypar_hg.num_hyperedges);
}

TEST_F(AKahyparPartitioner, ExtractsBlock1WithCorrectNumberOfHyperedges) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 1, node_mapping);
  ASSERT_EQ(3, kahypar_hg.num_hyperedges);
}

TEST_F(AKahyparPartitioner, ExtractsBlock0WithCorrectHyperedgeIndices) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 0, node_mapping);
  ASSERT_TRUE(isEqual(kahypar_hg.hyperedge_indices, {0, 2, 4, 5}));
}

TEST_F(AKahyparPartitioner, ExtractsBlock1WithCorrectHyperedgeIndices) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 1, node_mapping);
  ASSERT_TRUE(isEqual(kahypar_hg.hyperedge_indices, {0, 2, 5, 7}));
}

TEST_F(AKahyparPartitioner, ExtractsBlock0WithCorrectHyperedges) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 0, node_mapping);
  ASSERT_TRUE(isEqual(kahypar_hg.hyperedges, {0, 2, 0, 1, 2}));
}

TEST_F(AKahyparPartitioner, ExtractsBlock1WithCorrectHyperedges) {
  assignPartitionIDs(hypergraph);
  KaHyParHypergraph kahypar_hg = extractBlockAsKaHyParHypergraph(hypergraph, 1, node_mapping);
  ASSERT_TRUE(isEqual(kahypar_hg.hyperedges, {0, 1, 0, 1, 3, 2, 3}));
}

} // namespace mt_kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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



#include "gmock/gmock.h"
#include "tbb/task_group.h"

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/partition/fixed_vertices/fixed_vertex_removal.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {


class AFixedVertexRemoval : public Test {

 public:
  using Hypergraph = StaticHypergraph;
  using Factory = typename Hypergraph::Factory;

  AFixedVertexRemoval() :
    hypergraph(Factory::construct(
      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })) { }

  void addFixedVertices(const PartitionID k,
                        const vec<PartitionID>& fixed_vertex_blocks) {
    FixedVertexSupport<Hypergraph> fixed_vertices(hypergraph.initialNumNodes(), k);
    fixed_vertices.setHypergraph(&hypergraph);
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      if ( fixed_vertex_blocks[hn] != kInvalidPartition ) {
        fixed_vertices.fixToBlock(hn, fixed_vertex_blocks[hn]);
      }
    }
    hypergraph.addFixedVertexSupport(std::move(fixed_vertices));
  }

  void verifyFixedVertexFreeSubhypergraph(const Hypergraph& hypergraph,
                                          const vec<HypernodeID>& hn_mapping,
                                          const HypernodeID expected_num_nodes,
                                          const vec<vec<HypernodeID>>& expected_pins) {
    const HyperedgeID expected_num_edges = expected_pins.size();
    ASSERT_EQ(expected_num_nodes, hypergraph.initialNumNodes());
    ASSERT_EQ(expected_num_edges, hypergraph.initialNumEdges());

    vec<HypernodeID> to_original(expected_num_nodes, 0);
    for ( size_t i = 0; i < hn_mapping.size(); ++i ) {
      if ( hn_mapping[i] != kInvalidHypernode ) {
        to_original[hn_mapping[i]] = i;
      }
    }

    for ( const HyperedgeID& he : hypergraph.edges() ) {
      ASSERT_EQ(UL(hypergraph.edgeSize(he)), expected_pins[he].size());
      size_t i = 0;
      for ( const HypernodeID& pin : hypergraph.pins(he) ) {
        ASSERT_EQ(expected_pins[he][i++], to_original[pin]);
      }
    }
  }

  Hypergraph hypergraph;
};


TEST_F(AFixedVertexRemoval, ConstructFixedVertexFreeSubhypergraph1) {
  addFixedVertices(3, { kInvalidPartition, kInvalidPartition, 0, kInvalidPartition, 1, kInvalidPartition, 2 });
  ExtractedHypergraph<StaticHypergraph> extracted_hg = FixedVertexRemoval<StaticHypergraph>::remove(hypergraph);
  verifyFixedVertexFreeSubhypergraph(extracted_hg.hg, extracted_hg.hn_mapping, 4, { { 0, 1, 3 } });
}

TEST_F(AFixedVertexRemoval, ConstructFixedVertexFreeSubhypergraph2) {
  addFixedVertices(3, { 0, kInvalidPartition, kInvalidPartition, kInvalidPartition, kInvalidPartition, kInvalidPartition, 2 });
  ExtractedHypergraph<StaticHypergraph> extracted_hg = FixedVertexRemoval<StaticHypergraph>::remove(hypergraph);
  verifyFixedVertexFreeSubhypergraph(extracted_hg.hg, extracted_hg.hn_mapping,
    5, { { 1, 3, 4 }, { 3, 4 }, { 2, 5 } });
}

TEST_F(AFixedVertexRemoval, ConstructFixedVertexFreeSubhypergraph3) {
  addFixedVertices(3, { 0, 1, 2, kInvalidPartition, kInvalidPartition, 2, 2 });
  ExtractedHypergraph<StaticHypergraph> extracted_hg = FixedVertexRemoval<StaticHypergraph>::remove(hypergraph);
  verifyFixedVertexFreeSubhypergraph(extracted_hg.hg, extracted_hg.hn_mapping,
    2, { { 3, 4 }, { 3, 4 } });
}

TEST_F(AFixedVertexRemoval, ConstructFixedVertexFreeSubhypergraph4) {
  addFixedVertices(3, { kInvalidPartition, 1, kInvalidPartition, 2, 1, kInvalidPartition, 2 });
  ExtractedHypergraph<StaticHypergraph> extracted_hg = FixedVertexRemoval<StaticHypergraph>::remove(hypergraph);
  verifyFixedVertexFreeSubhypergraph(extracted_hg.hg, extracted_hg.hn_mapping,
    3, { { 0, 2 }, { 2, 5 } });
}

TEST_F(AFixedVertexRemoval, ConstructFixedVertexFreeSubhypergraph5) {
  addFixedVertices(3, { kInvalidPartition, 1, 0, 2, 1, kInvalidPartition, 2 });
  ExtractedHypergraph<StaticHypergraph> extracted_hg = FixedVertexRemoval<StaticHypergraph>::remove(hypergraph);
  verifyFixedVertexFreeSubhypergraph(extracted_hg.hg, extracted_hg.hn_mapping, 2, { });
}

}  // namespace ds
}  // namespace mt_kahypar

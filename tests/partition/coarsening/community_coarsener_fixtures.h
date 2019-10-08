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

#include "mt-kahypar/definitions.h"
#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;
using ::testing::Eq;
using ::testing::Le;

namespace mt_kahypar {

class ACommunityCoarsener : public ds::AHypergraph<2> {

 private:
  using Base = ds::AHypergraph<2>;

 public:
  using Base::TestStreamingHypergraph;
  using Base::TestHypergraph;

  ACommunityCoarsener() :
    Base(),
    hypergraph(),
    context() {
    hypergraph = construct_hypergraph( 16,
      { {0, 1}, {0, 1, 3}, {1, 2, 3}, {2, 3, 4}, {2, 4},
        {4, 5}, {4, 5, 7}, {5, 6, 7}, {6, 7, 8}, {6, 8},
        {8, 9}, {8, 9, 11}, {9, 10, 11}, {10, 11, 12}, {10, 12},
        {12, 13}, {12, 13, 15}, {13, 14, 15} },
      { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 },
      { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 } );
    hypergraph.setCommunityNodeMapping( { 0, 0, 1, 1 } );

    context.coarsening.max_allowed_node_weight = std::numeric_limits<HypernodeWeight>::max();
    context.coarsening.contraction_limit = 8;
  }

  TestHypergraph hypergraph;
  Context context;
};

template< class HyperGraph >
void assignPartitionIDs(HyperGraph& hypergraph) {
  using StreamingHyperGraph = typename HyperGraph::StreamingHypergraph;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    PartitionID part_id = StreamingHyperGraph::get_numa_node_of_vertex(hn);
    hypergraph.setPartInfo(hn, part_id);
  }
  hypergraph.updateGlobalPartInfos();
}

template< class HyperGraph >
HypernodeID currentNumNodes( HyperGraph& hypergraph  ) {
  HypernodeID num_nodes = 0;
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    unused(hn);
    ++num_nodes;
  }
  return num_nodes;
}

template< class HyperGraph >
HyperedgeID currentNumEdges( HyperGraph& hypergraph  ) {
  HyperedgeID num_edges = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    unused(he);
    ++num_edges;
  }
  return num_edges;
}

template< class HyperGraph >
HypernodeID currentNumPins( HyperGraph& hypergraph  ) {
  HypernodeID num_pins = 0;
  for ( const HypernodeID& he : hypergraph.edges() ) {
    num_pins += hypergraph.edgeSize(he);
  }
  return num_pins;
}

template <class Coarsener>
void doCoarsening(Coarsener& coarsener) {
  coarsener.disableRandomization();
  coarsener.coarsen();
}

template <class Coarsener, class Hypergraph>
void decreasesNumberOfPins(Coarsener& coarsener,
                           Hypergraph& hypergraph,
                           const size_t number_of_pins) {
  doCoarsening(coarsener);
  ASSERT_THAT(currentNumPins(hypergraph), Eq(number_of_pins));
}

template <class Coarsener, class HyperGraph>
void decreasesNumberOfHyperedges(Coarsener& coarsener,
                                 HyperGraph& hypergraph,
                                 const HyperedgeID num_hyperedges) {
  doCoarsening(coarsener);
  ASSERT_THAT(currentNumEdges(hypergraph), Eq(num_hyperedges));
}

template <class Coarsener, class HyperGraph>
void removesHyperedgesOfSizeOneDuringCoarsening(Coarsener& coarsener,
                                                HyperGraph& hypergraph,
                                                const std::vector<HyperedgeID>& single_node_hes) {
  doCoarsening(coarsener);
  for ( const HyperedgeID& he : single_node_hes ) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void reAddsHyperedgesOfSizeOneDuringUncoarsening(Coarsener& coarsener,
                                                 HyperGraph& hypergraph,
                                                 const std::vector<HyperedgeID>& single_node_hes) {
  doCoarsening(coarsener);
  for ( const HyperedgeID& he : single_node_hes ) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
  assignPartitionIDs(hypergraph);
  coarsener.uncoarsen();
  for ( const HyperedgeID& he : single_node_hes ) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void removesParallelHyperedgesDuringCoarsening(Coarsener& coarsener,
                                               HyperGraph& hypergraph,
                                               const std::vector<HyperedgeID>& parallel_hes) {
  doCoarsening(coarsener);
  for ( const HyperedgeID& he : parallel_hes ) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void updatesEdgeWeightOfRepresentativeHyperedgeOnParallelHyperedgeRemoval(Coarsener& coarsener,
                                                                          HyperGraph& hypergraph,
                                                                          const std::vector<std::pair<HyperedgeID, HyperedgeWeight>>& he_weights) {
  doCoarsening(coarsener);
  for ( const auto& he_weight : he_weights ) {
    HyperedgeID he = he_weight.first;
    HyperedgeWeight weight = he_weight.second;
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
    ASSERT_THAT(hypergraph.edgeWeight(he), Eq(weight));
  }
}

template <class Coarsener, class HyperGraph>
void restoresParallelHyperedgesDuringUncoarsening(Coarsener& coarsener,
                                                  HyperGraph& hypergraph,
                                                  const std::vector<HyperedgeID>& parallel_hes) {
  doCoarsening(coarsener);
  for ( const HyperedgeID& he : parallel_hes ) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
  assignPartitionIDs(hypergraph);
  coarsener.uncoarsen();
  for ( const HyperedgeID& he : parallel_hes ) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void doesNotCoarsenUntilCoarseningLimit(Coarsener& coarsener,
                                        HyperGraph& hypergraph,
                                        Context& context,
                                        const HypernodeID contraction_limit,
                                        const HypernodeWeight max_allowed_node_weight,
                                        const size_t expected_num_nodes) {
  context.coarsening.contraction_limit = contraction_limit;
  context.coarsening.max_allowed_node_weight = max_allowed_node_weight;
  doCoarsening(coarsener);
  for (const HypernodeID& hn : hypergraph.nodes()) {
    ASSERT_THAT(hypergraph.nodeWeight(hn), Le(context.coarsening.max_allowed_node_weight));
  }
  ASSERT_THAT(currentNumNodes(hypergraph), Eq(expected_num_nodes));
}

} // namespace mt_kahypar

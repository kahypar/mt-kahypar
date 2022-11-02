/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/multilevel_coarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_coarsener.h"

using ::testing::Test;
using ::testing::Eq;
using ::testing::Le;

namespace mt_kahypar {

using PartitionedHyperGraph = mt_kahypar::PartitionedHypergraph;

class ACoarsener : public Test {
 private:

 public:
  ACoarsener() :
    hypergraph(HypergraphFactory::construct(16, 18, { { 0, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 3, 4 }, { 2, 4 },
                { 4, 5 }, { 4, 5, 7 }, { 5, 6, 7 }, { 6, 7, 8 }, { 6, 8 },
                { 8, 9 }, { 8, 9, 11 }, { 9, 10, 11 }, { 10, 11, 12 }, { 10, 12 },
                { 12, 13 }, { 12, 13, 15 }, { 13, 14, 15 } })),
    context(),
    nullptr_refiner(nullptr) {
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hypergraph.setCommunityID(hn, hn / 4);
    }

    context.partition.k = 2;
    context.partition.mode = Mode::direct;
    context.partition.objective = Objective::km1;
    context.coarsening.max_allowed_node_weight = std::numeric_limits<HypernodeWeight>::max();
    context.coarsening.contraction_limit = 8;
    context.coarsening.use_adaptive_max_allowed_node_weight = false;
    context.coarsening.minimum_shrink_factor = 1.0;
    context.coarsening.maximum_shrink_factor = 4.0;
    context.refinement.max_batch_size = 5;
    context.setupPartWeights(hypergraph.totalWeight());
  }

  Hypergraph hypergraph;
  Context context;
  std::unique_ptr<IRefiner> nullptr_refiner;
};

template <class Hypergraph>
void assignPartitionIDs(Hypergraph& hypergraph) {
  for (const HypernodeID& hn : hypergraph.nodes()) {
    PartitionID part_id = 0;
    hypergraph.setNodePart(hn, part_id);
  }
}

template <class Hypergraph>
HypernodeID currentNumNodes(Hypergraph& hypergraph) {
  HypernodeID num_nodes = 0;
  for (const HypernodeID& hn : hypergraph.nodes()) {
    unused(hn);
    ++num_nodes;
  }
  return num_nodes;
}

template <class Hypergraph>
HyperedgeID currentNumEdges(Hypergraph& hypergraph) {
  HyperedgeID num_edges = 0;
  for (const HyperedgeID& he : hypergraph.edges()) {
    unused(he);
    ++num_edges;
  }
  return num_edges;
}

template <class Hypergraph>
HypernodeID currentNumPins(Hypergraph& hypergraph) {
  HypernodeID num_pins = 0;
  for (const HypernodeID& he : hypergraph.edges()) {
    num_pins += hypergraph.edgeSize(he);
  }
  return num_pins;
}

template <class Coarsener>
void doCoarsening(Coarsener& coarsener) {
  coarsener.disableRandomization();
  coarsener.coarsen();
}

template <class Coarsener>
void decreasesNumberOfPins(Coarsener& coarsener,
                           const size_t number_of_pins) {
  doCoarsening(coarsener);
  ASSERT_THAT(currentNumPins(coarsener.coarsestHypergraph()), Eq(number_of_pins));
}

template <class Coarsener>
void decreasesNumberOfHyperedges(Coarsener& coarsener,
                                 const HyperedgeID num_hyperedges) {
  doCoarsening(coarsener);
  ASSERT_THAT(currentNumEdges(coarsener.coarsestHypergraph()), Eq(num_hyperedges));
}

template <class Coarsener, class Hypergraph>
void removesHyperedgesOfSizeOneDuringCoarsening(Coarsener& coarsener,
                                                Hypergraph& hypergraph,
                                                const std::vector<HyperedgeID>& single_node_hes) {
  doCoarsening(coarsener);
  for (const HyperedgeID& he : single_node_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
}

template <class Coarsener, class Hypergraph>
void removesParallelHyperedgesDuringCoarsening(Coarsener& coarsener,
                                               Hypergraph& hypergraph,
                                               const std::vector<HyperedgeID>& parallel_hes) {
  doCoarsening(coarsener);
  for (const HyperedgeID& he : parallel_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
}

template <class Coarsener, class Hypergraph>
void updatesEdgeWeightOfRepresentativeHyperedgeOnParallelHyperedgeRemoval(Coarsener& coarsener,
                                                                          Hypergraph& hypergraph,
                                                                          const std::vector<std::pair<HyperedgeID, HyperedgeWeight> >& he_weights) {
  doCoarsening(coarsener);
  for (const auto& he_weight : he_weights) {
    HyperedgeID he = he_weight.first;
    HyperedgeWeight weight = he_weight.second;
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
    ASSERT_THAT(hypergraph.edgeWeight(he), Eq(weight));
  }
}

template <class Coarsener, class Hypergraph>
void doesNotCoarsenUntilCoarseningLimit(Coarsener& coarsener,
                                        Hypergraph& hypergraph,
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
}  // namespace mt_kahypar

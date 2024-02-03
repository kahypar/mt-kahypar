#pragma once

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

template <typename TypeTraits>
class TriangleCounter {
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  TriangleCounter(const Context& context, const Hypergraph& hypergraph): 
    _context(context), _original_hypergraph(hypergraph), _original_hypergraph_weights({}){ }

  TriangleCounter(const TriangleCounter&) = delete;
  TriangleCounter& operator=(const TriangleCounter&) = delete;

  TriangleCounter(TriangleCounter&&) = delete;
  TriangleCounter& operator=(TriangleCounter&&) = delete;

  // ! Counts the number of triangles in the hypergraph and multiplies the edge weights 
  // ! with the triangle count on each hyperedge respectively.
  void countTrianglesAndReplaceEdgeWeights(Hypergraph& hypergraph) {
    _original_hypergraph_weights.resize(hypergraph.initialNumEdges());
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      _original_hypergraph_weights[he] = hypergraph.edgeWeight(he);
      hypergraph.multEdgeWeightWithTriangleCount(he);
    });
    _triangle_graph = hypergraph.copy();
  }

  // ! Restores all previously changed hyperedges.
  void replaceInitialWeights(Hypergraph& hypergraph) {
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      hypergraph.setEdgeWeight(he, _original_hypergraph_weights[he]);
    });
  }

 private:
  const Context& _context;
  const Hypergraph& _original_hypergraph;
  std::vector<HypernodeID> _original_hypergraph_weights;
  Hypergraph _triangle_graph;
};

}  // namespace mt_kahypar

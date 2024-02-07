#pragma once

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

template <typename TypeTraits>
class TriangleCounter {
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  TriangleCounter(const Context& context, Hypergraph& hypergraph): 
    _context(context), _original_hypergraph(hypergraph), _original_hypergraph_weights({}) {
      _triangle_graph = hypergraph.copy();
      _original_hypergraph_weights.resize(hypergraph.initialNumEdges());
      _triangle_graph.doParallelForAllEdges([&](const HyperedgeID& he) {
        _original_hypergraph_weights[he] = hypergraph.edgeWeight(he);
        _triangle_graph.multEdgeWeightWithTriangleCount(he);
      });
    }

  TriangleCounter(const TriangleCounter&) = delete;
  TriangleCounter& operator=(const TriangleCounter&) = delete;

  TriangleCounter(TriangleCounter&&) = delete;
  TriangleCounter& operator=(TriangleCounter&&) = delete;

  Hypergraph& getTriangleGraph() {
    return _triangle_graph;
  }

  // ! Counts the number of triangles in the hypergraph and multiplies the edge weights 
  // ! with the triangle count on each hyperedge respectively.
  void replaceEdgeWeights() {
    _original_hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      _original_hypergraph.setEdgeWeight(he, _triangle_graph.edgeWeight(he));
    });
  }

  // ! Restores all previously changed hyperedges.
  void replaceInitialWeights() {
    _original_hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      _original_hypergraph.setEdgeWeight(he, _original_hypergraph_weights[he]);
    });
  }

 private:
  const Context& _context;
  Hypergraph& _original_hypergraph;
  std::vector<HypernodeID> _original_hypergraph_weights;
  Hypergraph _triangle_graph;
};

}  // namespace mt_kahypar

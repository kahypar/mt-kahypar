#pragma once

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

template <typename TypeTraits>
class TriangleCounter {
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  TriangleCounter(const Context& context, const Hypergraph& hypergraph): 
    _context(context), _original_hypergraph(hypergraph) { }

  TriangleCounter(const TriangleCounter&) = delete;
  TriangleCounter& operator=(const TriangleCounter&) = delete;

  TriangleCounter(TriangleCounter&&) = delete;
  TriangleCounter& operator=(TriangleCounter&&) = delete;

  // ! Counts the number of triangles in the hypergraph and multiplies the edge weights 
  // ! with the triangle count on each hyperedge respectively.
  void countTrianglesAndReplaceEdgeWeights(Hypergraph& hg) {
    for (const HyperedgeID& he : _original_hypergraph.edges()) {
      hg.multEdgeWeightWithTriangleCount(he);
    }
  }

  // ! Restores all previously changed hyperedges.
  void replaceInitialWeights(PartitionedHypergraph& hypergraph) {
    for (const HyperedgeID& he : hypergraph.edges()) {
      hypergraph.setEdgeWeight(he, _original_hypergraph.edgeWeight(he));
    }
  }

 private:
  const Context& _context;
  const Hypergraph& _original_hypergraph;
};

}  // namespace mt_kahypar

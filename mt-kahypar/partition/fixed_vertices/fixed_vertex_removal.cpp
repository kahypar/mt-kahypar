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

#include "mt-kahypar/partition/fixed_vertices/fixed_vertex_removal.h"

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

namespace mt_kahypar {

template<typename Hypergraph>
ExtractedHypergraph<Hypergraph> FixedVertexRemoval<Hypergraph>::remove(const Hypergraph& hypergraph) {
  ExtractedHypergraph<Hypergraph> extracted_hg;
  vec<HypernodeID>& hn_mapping = extracted_hg.hn_mapping;
  vec<HyperedgeID> he_mapping;

  parallel::TBBPrefixSum<HypernodeID> hn_mapping_prefix_sum(hn_mapping);
  parallel::TBBPrefixSum<HyperedgeID> he_mapping_prefix_sum(he_mapping);
  tbb::parallel_invoke([&] {
    hn_mapping.assign(hypergraph.initialNumNodes() + 1, 0);
    hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      if ( !hypergraph.isFixed(hn) ) {
        // Mark free vertices as valid
        hn_mapping[hn + 1] = 1;
      }
    });
    // Prefix sums computes mapping of the nodes of the input hypergraph to the nodes
    // of the fixed vertex free subhypergraph
    tbb::parallel_scan(tbb::blocked_range<size_t>(UL(0), UL(hypergraph.initialNumNodes() + 1)), hn_mapping_prefix_sum);
  }, [&] {
    he_mapping.assign(hypergraph.initialNumEdges() + 1, 0);
    hypergraph.doParallelForAllEdges([&](const HypernodeID he) {
      size_t edge_size = 0;
      for ( const HypernodeID& pin : hypergraph.pins(he) ) {
        if ( !hypergraph.isFixed(pin) ) {
          ++edge_size;
        }
      }
      if ( edge_size > 1 ) {
        // Mark edges with more than one pin as valid
        he_mapping[he + 1] = 1;
      }
    });
    // Prefix sums computes mapping of the edges of the input hypergraph to the edges
    // of the fixed vertex free subhypergraph
    tbb::parallel_scan(tbb::blocked_range<size_t>(UL(0), UL(hypergraph.initialNumEdges() + 1)), he_mapping_prefix_sum);
  });
  // Remove sentinel
  hn_mapping.pop_back();

  const HypernodeID num_nodes = hn_mapping_prefix_sum.total_sum();
  const HyperedgeID num_edges = he_mapping_prefix_sum.total_sum();
  HyperedgeVector edge_vector;
  vec<HyperedgeWeight> hyperedge_weights;
  vec<HypernodeWeight> hypernode_weights;
  tbb::parallel_invoke([&] {
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_edges);
    }, [&] {
      hyperedge_weights.resize(num_edges);
    });

    // Construct edge vector of the fixed vertex free subhypergraph
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      if ( he_mapping_prefix_sum.value(he + 1) ) {
        // Hyperedge contains more than one pin in fixed vertex free subhypergraph
        const HyperedgeID e = he_mapping[he];
        hyperedge_weights[e] = hypergraph.edgeWeight(he);
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          if ( !hypergraph.isFixed(pin) ) {
            edge_vector[e].push_back(hn_mapping[pin]);
          }
        }
      }
    });
  }, [&] {
    hypernode_weights.resize(num_nodes);
    hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( !hypergraph.isFixed(hn) ) {
        hypernode_weights[hn_mapping[hn]] = hypergraph.nodeWeight(hn);
      } else {
        hn_mapping[hn] = kInvalidHypernode;
      }
    });
  });

  // Construct fixed vertex free subhypergraph
  extracted_hg.hg = Factory::construct(num_nodes, num_edges,
    edge_vector, hyperedge_weights.data(), hypernode_weights.data(), true);
  return extracted_hg;
}

INSTANTIATE_CLASS_WITH_HYPERGRAPHS(FixedVertexRemoval)
}

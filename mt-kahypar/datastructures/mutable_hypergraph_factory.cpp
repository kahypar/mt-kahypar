/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mutable_hypergraph_factory.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::ds {

  MutableHypergraph MutableHypergraphFactory::construct(
          const HypernodeID num_hypernodes,
          const HyperedgeID num_hyperedges,
          const HyperedgeVector& edge_vector,
          const HyperedgeWeight* hyperedge_weight,
          const HypernodeWeight* hypernode_weight,
          [[maybe_unused]] const bool stable_construction_of_incident_edges) {
    MutableHypergraph hypergraph;
    hypergraph._num_hypernodes = num_hypernodes;
    hypergraph._num_hyperedges = num_hyperedges;
    hypergraph._hypernodes.resize(num_hypernodes);
    hypergraph._hyperedges.resize(num_hyperedges);

    hypergraph._num_pins = 0;
    hypergraph._total_degree = 0;
    hypergraph._max_edge_size = 0;

    //TODO is this equivalent to stable_construction ?


    // Create Hypernode adjacency structure
    std::vector<std::vector<HyperedgeID>> incident_nets(num_hypernodes);

    for (HyperedgeID i = 0; i < num_hyperedges; ++i) {
      hypergraph._hyperedges[i] = MutableHypergraph::Hyperedge(edge_vector[i],
                                                               hyperedge_weight ? hyperedge_weight[i] : 1);
      for (const HypernodeID& pin : edge_vector[i]) {
        ASSERT(pin < num_hypernodes, V(pin) << V(num_hypernodes));
        incident_nets[pin].push_back(i);
      }

      hypergraph._num_pins += edge_vector[i].size();
      if ( hypergraph._max_edge_size < edge_vector[i].size() ) {
        hypergraph._max_edge_size = edge_vector[i].size();
      }
    }

    for (HypernodeID i = 0; i < num_hypernodes; ++i) {
      hypergraph._hypernodes[i] = MutableHypergraph::Hypernode(incident_nets[i], hypernode_weight ? hypernode_weight[i] : 1);
      hypergraph._total_degree += hypergraph._hypernodes[i].size();
      //TODO total weight?
    }

    hypergraph._community_ids.resize(num_hypernodes, 0);

    hypergraph.computeAndSetTotalNodeWeight(parallel_tag_t());

    return hypergraph;
  }

}
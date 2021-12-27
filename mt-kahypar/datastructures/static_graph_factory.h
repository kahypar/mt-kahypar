/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#pragma once

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/datastructures/static_graph.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"


namespace mt_kahypar {
namespace ds {

class StaticGraphFactory {
  using EdgeVector = parallel::scalable_vector<std::pair<HypernodeID, HypernodeID>>;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using Counter = parallel::scalable_vector<size_t>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;

 public:
  static StaticGraph construct(const HypernodeID num_nodes,
                               const HyperedgeID num_edges,
                               const HyperedgeVector& edge_vector,
                               const HyperedgeWeight* edge_weight = nullptr,
                               const HypernodeWeight* node_weight = nullptr,
                               const bool stable_construction_of_incident_edges = false);

  // ! Provides a more performant construction method by using continuous space for the edges
  // ! (instead of a separate vec per edge).
  // ! No backwards edges allowed, i.e. each edge is unique
  static StaticGraph construct_from_graph_edges(const HypernodeID num_nodes,
                                                const HyperedgeID num_edges,
                                                const EdgeVector& edge_vector,
                                                const HyperedgeWeight* edge_weight = nullptr,
                                                const HypernodeWeight* node_weight = nullptr,
                                                const bool stable_construction_of_incident_edges = false);

  static std::pair<StaticGraph, parallel::scalable_vector<HypernodeID> > compactify(const StaticGraph&) {
    ERROR("Compactify not implemented for static graph.");
  }

 private:
  StaticGraphFactory() { }

  static void sort_incident_edges(StaticGraph& graph);
};

} // namespace ds
} // namespace mt_kahypar
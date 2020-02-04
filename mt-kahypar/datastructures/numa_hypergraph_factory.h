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

#pragma once

#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/numa_hypergraph.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template <typename Hypergraph = Mandatory,
          typename Factory = Mandatory,
          typename HardwareTopology = Mandatory,
          typename TBBNumaArena = Mandatory>
class NumaHypergraphFactory {

  using NumaHyperGraph = NumaHypergraph<Hypergraph, HardwareTopology, TBBNumaArena>;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;

 public:
  static NumaHyperGraph construct(const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HyperedgeVector& edge_vector,
                                    const HyperedgeWeight* hyperedge_weight = nullptr,
                                    const HypernodeWeight* hypernode_weight = nullptr) {
    NumaHyperGraph hypergraph;

    return hypergraph;
  }

 private:
  NumaHypergraphFactory() { }
};

} // namespace ds
} // namespace mt_kahypar
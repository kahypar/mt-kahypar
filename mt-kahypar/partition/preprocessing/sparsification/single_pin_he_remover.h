/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template<typename TypeTraits>
class SinglePinHyperedgeRemoverT {

  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::PartitionedHyperGraph;

 public:
  SinglePinHyperedgeRemoverT() :
    _removed_hes() { }

  SinglePinHyperedgeRemoverT(const SinglePinHyperedgeRemoverT&) = delete;
  SinglePinHyperedgeRemoverT & operator= (const SinglePinHyperedgeRemoverT &) = delete;

  SinglePinHyperedgeRemoverT(SinglePinHyperedgeRemoverT&&) = delete;
  SinglePinHyperedgeRemoverT & operator= (SinglePinHyperedgeRemoverT &&) = delete;

  // ! Removes single-pin hyperedges since they do not contribute
  // ! to cut or km1 metric.
  HyperedgeID removeSinglePinHyperedges(HyperGraph& hypergraph) {
    HyperedgeID num_removed_single_node_hes = 0;
    for (const HyperedgeID& he : hypergraph.edges()) {
      if (hypergraph.edgeSize(he) == 1) {
        ++num_removed_single_node_hes;
        hypergraph.removeEdge(he);
        _removed_hes.push_back(he);
      }
    }
    return num_removed_single_node_hes;
  }

  // ! Restore single-pin hyperedges
  void restoreSinglePinHyperedges(PartitionedHyperGraph& hypergraph) {
    for (const HyperedgeID& he : _removed_hes) {
      hypergraph.restoreSinglePinHyperedge(he);
    }
  }

 private:
  parallel::scalable_vector<HyperedgeID> _removed_hes;
};

using SinglePinHyperedgeRemover = SinglePinHyperedgeRemoverT<GlobalTypeTraits>;

}  // namespace mt_kahypar

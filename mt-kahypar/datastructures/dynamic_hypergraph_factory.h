/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/enumerable_thread_specific.h"


#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace ds {

class DynamicHypergraphFactory {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using Counter = parallel::scalable_vector<size_t>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;
  using ThreadLocalBitset = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<>>;
  using ThreadLocalBitvector = tbb::enumerable_thread_specific<parallel::scalable_vector<bool>>;


 public:
  static DynamicHypergraph construct(const TaskGroupID task_group_id,
                                    const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HyperedgeVector& edge_vector,
                                    const HyperedgeWeight* hyperedge_weight = nullptr,
                                    const HypernodeWeight* hypernode_weight = nullptr,
                                    const bool stable_construction_of_incident_edges = false);

  /**
   * Compactifies a given hypergraph such that it only contains enabled vertices and hyperedges within
   * a consecutive range of IDs.
   */
  static std::pair<DynamicHypergraph,
                   parallel::scalable_vector<HypernodeID> > compactify(const TaskGroupID task_group_id,
                                                                       const DynamicHypergraph& hypergraph);

 private:
  DynamicHypergraphFactory() { }
};

} // namespace ds
} // namespace mt_kahypar
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

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"


namespace mt_kahypar {
namespace ds {

class StaticHypergraphFactory {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using Counter = parallel::scalable_vector<size_t>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;

 public:
  static StaticHypergraph construct(const TaskGroupID task_group_id,
                                    const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HyperedgeVector& edge_vector,
                                    const HyperedgeWeight* hyperedge_weight = nullptr,
                                    const HypernodeWeight* hypernode_weight = nullptr,
                                    const bool stable_construction_of_incident_edges = false);

  static std::pair<StaticHypergraph,
          parallel::scalable_vector<HypernodeID> > compactify(const TaskGroupID ,
                                                              const StaticHypergraph&)
                                                              { ERROR("Compactify not implemented for static hypergraph."); }

 private:
  StaticHypergraphFactory() { }
};

} // namespace ds
} // namespace mt_kahypar
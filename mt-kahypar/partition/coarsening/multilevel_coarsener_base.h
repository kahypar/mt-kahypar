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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/coarsening/level.h"


namespace mt_kahypar {

class MultilevelCoarsenerBase {
 private:

  static constexpr bool debug = false;

 public:
  MultilevelCoarsenerBase(Hypergraph& hypergraph,
                          const Context& context,
                          const bool top_level) :
          _is_finalized(false),
          _hg(hypergraph),
          _context(context),
          _top_level(top_level) {}

  MultilevelCoarsenerBase(const MultilevelCoarsenerBase&) = delete;
  MultilevelCoarsenerBase(MultilevelCoarsenerBase&&) = delete;
  MultilevelCoarsenerBase & operator= (const MultilevelCoarsenerBase &) = delete;
  MultilevelCoarsenerBase & operator= (MultilevelCoarsenerBase &&) = delete;

  virtual ~MultilevelCoarsenerBase() noexcept {
/*
    tbb::parallel_for(0UL, _hierarchy->size(), [&](const size_t i) {
      (*_hierarchy)[i].freeInternalData();
    }, tbb::static_partitioner());
*/
  }

 protected:

  HypernodeID currentNumNodes() const {
    if ( _hierarchy->empty() ) {
      return _hg.initialNumNodes();
    } else {
      return _hierarchy->back().contractedHypergraph().initialNumNodes();
    }
  }

  Hypergraph& currentHypergraph() {
    if ( _hierarchy->empty() ) {
      return _hg;
    } else {
      return _hierarchy->back().contractedHypergraph();
    }
  }

  PartitionedHypergraph& currentPartitionedHypergraph() {
    ASSERT(_is_finalized);
    if ( _hierarchy->empty() ) {
      return *_partitioned_hg;
    } else {
      return _hierarchy->back().contractedPartitionedHypergraph();
    }
  }

  void finalize();

  void performMultilevelContraction(
          parallel::scalable_vector<HypernodeID>&& communities,
          const HighResClockTimepoint& round_start);

 protected:
  bool _is_finalized;
  Hypergraph& _hg;
  std::shared_ptr<PartitionedHypergraph> _partitioned_hg;
  const Context& _context;
  const bool _top_level;
  std::shared_ptr<vec<Level>> _hierarchy;
};
}  // namespace mt_kahypar

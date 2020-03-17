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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {

namespace {
  using Hyperedge = parallel::scalable_vector<HypernodeID>;
} // namespace

class UnionCombiner {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline Hyperedge combine(const Hyperedge& lhs, const Hyperedge& rhs) {
    ASSERT(std::is_sorted(lhs.begin(), lhs.end()));
    ASSERT(std::is_sorted(rhs.begin(), rhs.end()));
    Hyperedge combined_he;
    size_t i = 0;
    size_t j = 0;
    while ( i < lhs.size() && j < rhs.size() ) {
      if ( lhs[i] == rhs[j] ) {
        combined_he.push_back(lhs[i]);
        ++i;
        ++j;
      } else if ( lhs[i] < rhs[j] ) {
        combined_he.push_back(lhs[i]);
        ++i;
      } else {
        combined_he.push_back(rhs[j]);
        ++j;
      }
    }

    ASSERT(i == lhs.size() || j == rhs.size());
    for ( ; i < lhs.size(); ++i ) {
      combined_he.push_back(lhs[i]);
    }
    for ( ; j < rhs.size(); ++j ) {
      combined_he.push_back(rhs[j]);
    }
    ASSERT(std::is_sorted(combined_he.begin(), combined_he.end()));
    return combined_he;
  }
};

class SamplingCombiner {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline Hyperedge combine(const Hyperedge& lhs, const Hyperedge& rhs) {
    Hyperedge combined_he;
    int cpu_id = sched_getcpu();
    for ( const HypernodeID& pin : lhs ) {
      if ( utils::Randomize::instance().flipCoin(cpu_id) ) {
        combined_he.push_back(pin);
      }
    }
    for ( const HypernodeID& pin : rhs ) {
      if ( utils::Randomize::instance().flipCoin(cpu_id) ) {
        combined_he.push_back(pin);
      }
    }
    std::sort(combined_he.begin(), combined_he.end());
    combined_he.erase(std::unique(combined_he.begin(), combined_he.end()), combined_he.end());
    return combined_he;
  }
};

class UndefinedCombiner {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline Hyperedge combine(const Hyperedge& lhs, const Hyperedge&) {
    Hyperedge combined_he(lhs);
    ERROR("Similiar net combine strategy is undefined");
    return combined_he;
  }
};

}  // namespace mt_kahypar

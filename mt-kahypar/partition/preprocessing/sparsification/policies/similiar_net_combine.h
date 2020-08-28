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
  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Hyperedge combine(
    const Hypergraph&, const Hyperedge& lhs, const Hyperedge& rhs) {
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

class MaxSizeCombiner {
 public:
  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Hyperedge combine(
    const Hypergraph&, const Hyperedge& lhs, const Hyperedge& rhs) {
    if ( lhs.size() < rhs.size() ) {
      return rhs;
    } else {
      return lhs;
    }
  }
};

class NetImportanceCombiner {
 public:
  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Hyperedge combine(
    const Hypergraph& hypergraph, const Hyperedge& lhs, const Hyperedge& rhs) {
    HyperedgeWeight score_lhs = 0;
    for ( const HypernodeID& pin : lhs ) {
      for ( const HyperedgeID& he : hypergraph.incidentEdges(pin) ) {
        score_lhs += hypergraph.edgeWeight(he);
      }
    }

    HyperedgeWeight score_rhs = 0;
    for ( const HypernodeID& pin : rhs ) {
      for ( const HyperedgeID& he : hypergraph.incidentEdges(pin) ) {
        score_rhs += hypergraph.edgeWeight(he);
      }
    }

    if ( score_lhs > score_rhs ) {
      return lhs;
    } else {
      return rhs;
    }
  }
};


class UndefinedCombiner {
 public:
  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Hyperedge combine(
    const Hypergraph&, const Hyperedge& lhs, const Hyperedge&) {
    Hyperedge combined_he(lhs);
    ERROR("Similiar net combine strategy is undefined");
    return combined_he;
  }
};

}  // namespace mt_kahypar

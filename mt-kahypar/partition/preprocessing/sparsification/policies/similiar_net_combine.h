/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

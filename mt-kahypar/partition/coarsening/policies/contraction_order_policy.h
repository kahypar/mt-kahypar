/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

using ContractionOrderPair = std::pair<HypernodeID, HypernodeID>;

class RatingOrder final : public kahypar::meta::PolicyBase {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline ContractionOrderPair order(const Hypergraph&,
                                                                           const HypernodeID u,
                                                                           const HypernodeID v) {
    return std::make_pair(u, v);
  }
};

class DegreeOrder final : public kahypar::meta::PolicyBase {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline ContractionOrderPair order(const Hypergraph& hypergraph,
                                                                           const HypernodeID u,
                                                                           const HypernodeID v) {
    if ( hypergraph.nodeDegree(u) < hypergraph.nodeDegree(v) ) {
      return std::make_pair(v, u);
    } else {
      return std::make_pair(u, v);
    }
  }
};

class HeavyNodeFirst final : public kahypar::meta::PolicyBase {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline ContractionOrderPair order(const Hypergraph& hypergraph,
                                                                           const HypernodeID u,
                                                                           const HypernodeID v) {
    if ( hypergraph.nodeWeight(u) < hypergraph.nodeWeight(v) ) {
      return std::make_pair(v, u);
    } else {
      return std::make_pair(u, v);
    }
  }
};

class HeavyNodeSecond final : public kahypar::meta::PolicyBase {
 public:
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline ContractionOrderPair order(const Hypergraph& hypergraph,
                                                                           const HypernodeID u,
                                                                           const HypernodeID v) {
    if ( hypergraph.nodeWeight(u) < hypergraph.nodeWeight(v) ) {
      return std::make_pair(u, v);
    } else {
      return std::make_pair(v, u);
    }
  }
};

using ContractionOrderPolicies = kahypar::meta::Typelist<RatingOrder, DegreeOrder, HeavyNodeFirst, HeavyNodeSecond>;
}  // namespace mt_kahypar

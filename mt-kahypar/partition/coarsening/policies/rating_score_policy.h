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
class HeavyEdgeScore final : public kahypar::meta::PolicyBase {
 public:
  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline RatingType score(const HyperGraph& hypergraph,
                                                                 const HyperedgeID he,
                                                                 const PartitionID community_id) {
    return static_cast<RatingType>(hypergraph.edgeWeight(he, community_id)) / (hypergraph.edgeSize(he) - 1);
  }
};

using RatingScorePolicies = kahypar::meta::Typelist<HeavyEdgeScore>;
}  // namespace mt_kahypar

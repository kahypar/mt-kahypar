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

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
class VertexObjectivePolicy final : public kahypar::meta::PolicyBase {
 public:
  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline HypernodeID objective(const HyperGraph& hypergraph, const PartitionID community) {
    return hypergraph.numCommunityHypernodes(community);
  }

  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline HypernodeID total(const HyperGraph& hypergraph) {
    return hypergraph.initialNumNodes();
  }
};

class VertexDegreeObjectivePolicy final : public kahypar::meta::PolicyBase {
 public:
  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline HypernodeID objective(const HyperGraph& hypergraph, const PartitionID community) {
    return hypergraph.communityDegree(community);
  }

  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline HypernodeID total(const HyperGraph& hypergraph) {
    return hypergraph.initialNumPins();
  }
};

class PinObjectivePolicy final : public kahypar::meta::PolicyBase {
 public:
  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline HypernodeID objective(const HyperGraph& hypergraph, const PartitionID community) {
    return hypergraph.numCommunityPins(community);
  }

  template <typename HyperGraph>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static inline HypernodeID total(const HyperGraph& hypergraph) {
    return hypergraph.initialNumPins();
  }
};

using ObjectivePolicyClasses = kahypar::meta::Typelist<VertexObjectivePolicy, VertexDegreeObjectivePolicy, PinObjectivePolicy>;
}  // namespace mt_kahypar

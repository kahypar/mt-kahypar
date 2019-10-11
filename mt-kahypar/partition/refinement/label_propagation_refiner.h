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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

namespace mt_kahypar {

template< typename TypeTraits,
          typename ExecutionPolicy = Mandatory >
class LabelPropagationRefinerT final : public IRefiner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool debug = false;

 public:

  explicit LabelPropagationRefinerT(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context),
    _current_level(0),
    _execution_policy() {
    initialize();
  }

  LabelPropagationRefinerT(const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT(LabelPropagationRefinerT&&) = delete;

  LabelPropagationRefinerT& operator= (const LabelPropagationRefinerT&) = delete;
  LabelPropagationRefinerT& operator= (LabelPropagationRefinerT&&) = delete;

  ~LabelPropagationRefinerT() override = default;

 private:
  bool refineImpl(std::vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& best_metrics) override final {
    unused(refinement_nodes);
    unused(best_metrics);
    ++_current_level;
    if ( !_execution_policy.execute(_current_level) ) {
      return false;
    }

    DBG << V(_context.refinement.label_propagation.algorithm) << V(_current_level);
    return false;
  }

  void initialize() {
    HypernodeID current_num_nodes = 0;
    for ( const HypernodeID& hn : _hg.nodes() ) {
      unused(hn);
      ++current_num_nodes;
    }
    _execution_policy.initialize(_hg, current_num_nodes);
  }

  HyperGraph& _hg;
  const Context& _context;
  size_t _current_level;
  ExecutionPolicy _execution_policy;
};

template< typename ExecutionPolicy = Mandatory >
using LabelPropagationRefiner = LabelPropagationRefinerT<GlobalTypeTraits, ExecutionPolicy>;

}  // namespace kahypar

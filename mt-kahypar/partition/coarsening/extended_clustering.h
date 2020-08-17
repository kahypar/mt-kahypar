/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "multilevel_coarsener_base.h"
#include "i_coarsener.h"

namespace mt_kahypar {

class ExtendedClustering : public ICoarsener, private MultilevelCoarsenerBase {
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static constexpr HypernodeID invalid_node = std::numeric_limits<HypernodeID>::max();
  using Base = MultilevelCoarsenerBase;

public:

  ExtendedClustering(Hypergraph& hypergraph,
  const Context& context,
  const TaskGroupID task_group_id,
  const bool top_level) : Base(hypergraph, context, task_group_id, top_level) { }

private:

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
  }

  PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation, std::unique_ptr<IRefiner>& fm) override ;

  void coarsenImpl() override ;
};

/*
  void ExtendedClustering::coarsenImpl() {

  }
*/

}
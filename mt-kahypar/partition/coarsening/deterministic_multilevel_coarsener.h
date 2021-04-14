/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/utils/reproducible_random.h"

namespace mt_kahypar {
class DeterministicMultilevelCoarsener :  public ICoarsener,
                                          private MultilevelCoarsenerBase
{
public:
  DeterministicMultilevelCoarsener(Hypergraph& hypergraph, const Context& context, const TaskGroupID task_group_id,
                                   const bool top_level) :
    Base(hypergraph, context, task_group_id, top_level)
  {
  }

  ~DeterministicMultilevelCoarsener() {

  }

private:

  using Base = MultilevelCoarsenerBase;
  using Base::_context;
  using Base::_task_group_id;

  utils::ParallelPermutation<HypernodeID> permutation;

  void coarsenImpl() override;

  PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                        std::unique_ptr<IRefiner>& fm) override {
    return Base::doUncoarsen(label_propagation, fm);
  }

  Hypergraph& coarsestHypergraphImpl() override {
    return Base::currentHypergraph();
  }

  PartitionedHypergraph& coarsestPartitionedHypergraphImpl() override {
    return Base::currentPartitionedHypergraph();
  }

};
}
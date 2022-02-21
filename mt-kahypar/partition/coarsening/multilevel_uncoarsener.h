/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/partition/coarsening/i_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/uncoarsener_base.h"
namespace mt_kahypar {

  class MultilevelUncoarsener : public IUncoarsener,
                                private UncoarsenerBase {

  public:
    MultilevelUncoarsener(Hypergraph& hypergraph,
                        const Context& context,
                        UncoarseningData& uncoarseningData) :
      UncoarsenerBase(hypergraph, context, uncoarseningData) { }

  MultilevelUncoarsener(const MultilevelUncoarsener&) = delete;
  MultilevelUncoarsener(MultilevelUncoarsener&&) = delete;
  MultilevelUncoarsener & operator= (const MultilevelUncoarsener &) = delete;
  MultilevelUncoarsener & operator= (MultilevelUncoarsener &&) = delete;

  private:
  PartitionedHypergraph&& doUncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                      std::unique_ptr<IRefiner>& fm);

  void refine(
    PartitionedHypergraph& partitioned_hypergraph,
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm,
    std::unique_ptr<IRefiner>& flows,
    Metrics& current_metrics,
    const double time_limit);

  PartitionedHypergraph&& uncoarsenImpl(
    std::unique_ptr<IRefiner>& label_propagation,
    std::unique_ptr<IRefiner>& fm) override {
    return doUncoarsen(label_propagation, fm);
  }
  };

}

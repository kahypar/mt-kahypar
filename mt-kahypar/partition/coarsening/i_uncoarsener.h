/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
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

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

namespace mt_kahypar {

  class IUncoarsener {

  public:
    IUncoarsener(const IUncoarsener&) = delete;
    IUncoarsener(IUncoarsener&&) = delete;
    IUncoarsener & operator= (const IUncoarsener &) = delete;
    IUncoarsener & operator= (IUncoarsener &&) = delete;

    PartitionedHypergraph&& uncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                      std::unique_ptr<IRefiner>& fm) {
      return uncoarsenImpl(label_propagation, fm);
    }

    virtual ~IUncoarsener() = default;

  protected:
    IUncoarsener() = default;

  private:
    virtual PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                                  std::unique_ptr<IRefiner>& fm) = 0;
  };
}

/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"

namespace mt_kahypar {

class DeepInitialPartitioner: public IInitialPartitioner {
 private:
  static constexpr bool enable_heavy_assert = false;

 public:
  DeepInitialPartitioner(PartitionedHypergraph& hypergraph, const Context& context);
  DeepInitialPartitioner(const DeepInitialPartitioner&) = delete;
  DeepInitialPartitioner(DeepInitialPartitioner&&) = delete;
  DeepInitialPartitioner & operator= (const DeepInitialPartitioner &) = delete;
  DeepInitialPartitioner & operator= (DeepInitialPartitioner &&) = delete;

 private:
  void initialPartitionImpl() final;

 private:
  PartitionedHypergraph& _hg;
  const Context& _context;
};

}  // namespace mt_kahypar

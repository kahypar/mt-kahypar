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

#include <libkahypar.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/kahypar.h"

namespace mt_kahypar {

template< typename TypeTraits >
class RecursiveInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool debug = false;

 public:
  RecursiveInitialPartitionerT(HyperGraph& hypergraph, const Context& context) :
    _hg(hypergraph),
    _context(context) { }

  RecursiveInitialPartitionerT(const RecursiveInitialPartitionerT&) = delete;
  RecursiveInitialPartitionerT(RecursiveInitialPartitionerT&&) = delete;
  RecursiveInitialPartitionerT& operator= (const RecursiveInitialPartitionerT&) = delete;
  RecursiveInitialPartitionerT& operator= (RecursiveInitialPartitionerT&&) = delete;

 private:
  void initialPartitionImpl() override final {
    LOG << "HELLO WORLD!";
  }

 private:
  HyperGraph& _hg;
  const Context& _context;
};

using RecursiveInitialPartitioner = RecursiveInitialPartitionerT<GlobalTypeTraits>;

}  // namespace mt_kahypar

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

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"

namespace mt_kahypar {
namespace multilevel {

static inline void partition(Hypergraph& hypergraph, const Context& context) {

  // ################## COARSENING ##################
  mt_kahypar::io::printCoarseningBanner(context);
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  std::unique_ptr<ICoarsener> coarsener =
    CoarsenerFactory::getInstance().createObject(
      context.coarsening.algorithm, hypergraph, context);
  coarsener->coarsen();
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("coarsening", "Coarsening",
    "", mt_kahypar::utils::Timer::Type::COARSENING, 2, std::chrono::duration<double>(end - start).count());

  if ( context.partition.verbose_output ) {
    mt_kahypar::io::printHypergraphInfo(hypergraph, "Coarsened Hypergraph");
  }

  // ################## INITIAL PARTITIONING ##################
  io::printInitialPartitioningBanner(context);
  start = std::chrono::high_resolution_clock::now();
  std::unique_ptr<IInitialPartitioner> initial_partitioner =
    InitialPartitionerFactory::getInstance().createObject(
      context.initial_partitioning.mode, hypergraph, context);
  initial_partitioner->initialPartition();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initial_partitioning", "Initial Partitioning",
    "", mt_kahypar::utils::Timer::Type::INITIAL_PARTITIONING, 3, std::chrono::duration<double>(end - start).count(), true);

  io::printPartitioningResults(hypergraph, context, "Initial Partitioning Results:");

  // ################## LOCAL SEARCH ##################
  io::printLocalSearchBanner(context);
  start = std::chrono::high_resolution_clock::now();
  std::unique_ptr<IRefiner> label_propagation =
    LabelPropagationFactory::getInstance().createObject(
      context.refinement.label_propagation.algorithm, hypergraph, context);

  coarsener->uncoarsen(label_propagation);
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("refinement", "Refinement",
    "", mt_kahypar::utils::Timer::Type::REFINEMENT, 4, std::chrono::duration<double>(end - start).count(), true);

  io::printPartitioningResults(hypergraph, context, "Local Search Results:");
}
}  // namespace multilevel
}  // namespace mt_kahypar

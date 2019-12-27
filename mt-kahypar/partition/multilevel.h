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

#include <memory>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

namespace mt_kahypar {
namespace multilevel {
static inline void partition(Hypergraph& hypergraph,
                             const Context& context,
                             const bool top_level,
                             const TaskGroupID task_group_id) {
  // ################## COARSENING ##################
  mt_kahypar::io::printCoarseningBanner(context);
  utils::Timer::instance().start_timer("coarsening", "Coarsening");
  std::unique_ptr<ICoarsener> coarsener =
    CoarsenerFactory::getInstance().createObject(
      context.coarsening.algorithm, hypergraph, context, task_group_id);
  coarsener->coarsen();
  utils::Timer::instance().stop_timer("coarsening");

  if (context.partition.verbose_output) {
    mt_kahypar::io::printHypergraphInfo(hypergraph, "Coarsened Hypergraph");
  }

  // ################## INITIAL PARTITIONING ##################
  io::printInitialPartitioningBanner(context);
  utils::Timer::instance().start_timer("initial_partitioning", "Initial Partitioning");
  std::unique_ptr<IInitialPartitioner> initial_partitioner =
    InitialPartitionerFactory::getInstance().createObject(
      context.initial_partitioning.mode, hypergraph, context, top_level, task_group_id);
  initial_partitioner->initialPartition();
  utils::Timer::instance().stop_timer("initial_partitioning");

  io::printPartitioningResults(hypergraph, context, "Initial Partitioning Results:");

  // ################## LOCAL SEARCH ##################
  io::printLocalSearchBanner(context);
  utils::Timer::instance().start_timer("refinement", "Refinement");
  std::unique_ptr<IRefiner> label_propagation =
    LabelPropagationFactory::getInstance().createObject(
      context.refinement.label_propagation.algorithm, hypergraph, context, task_group_id);

  coarsener->uncoarsen(label_propagation);
  utils::Timer::instance().stop_timer("refinement");

  io::printPartitioningResults(hypergraph, context, "Local Search Results:");
}
}  // namespace multilevel
}  // namespace mt_kahypar

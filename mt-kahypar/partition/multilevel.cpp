/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/multilevel.h"

#include <memory>

#include "tbb/task.h"

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"
#include "mt-kahypar/partition/deep_multilevel.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/coarsening/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar::multilevel {

namespace {
  void disableTimerAndStats(const Context& context) {
    if ( context.type == ContextType::main && context.partition.mode == Mode::direct ) {
      utils::Utilities& utils = utils::Utilities::instance();
      parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
      utils.getTimer(context.utility_id).disable();
      utils.getStats(context.utility_id).disable();
    }
  }

  void enableTimerAndStats(const Context& context) {
    if ( context.type == ContextType::main && context.partition.mode == Mode::direct ) {
      utils::Utilities& utils = utils::Utilities::instance();
      parallel::MemoryPool::instance().activate_unused_memory_allocations();
      utils.getTimer(context.utility_id).enable();
      utils.getStats(context.utility_id).enable();
    }
  }

  PartitionedHypergraph multilevel_partitioning(Hypergraph& hypergraph,
                                                const Context& context,
                                                const bool is_vcycle) {
    PartitionedHypergraph partitioned_hg;

    // ################## COARSENING ##################
    mt_kahypar::io::printCoarseningBanner(context);

    const bool nlevel = context.coarsening.algorithm == CoarseningAlgorithm::nlevel_coarsener;
    UncoarseningData uncoarseningData(nlevel, hypergraph, context);

    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
    timer.start_timer("coarsening", "Coarsening");
    {
      std::unique_ptr<ICoarsener> coarsener = CoarsenerFactory::getInstance().createObject(
        context.coarsening.algorithm, hypergraph, context, uncoarseningData);
      coarsener->coarsen();

      if (context.partition.verbose_output) {
        Hypergraph& coarsestHypergraph = coarsener->coarsestHypergraph();
        mt_kahypar::io::printHypergraphInfo(coarsestHypergraph,
          "Coarsened Hypergraph", context.partition.show_memory_consumption);
      }
    }
    timer.stop_timer("coarsening");

    // ################## INITIAL PARTITIONING ##################
    io::printInitialPartitioningBanner(context);
    timer.start_timer("initial_partitioning", "Initial Partitioning");
    PartitionedHypergraph& phg = uncoarseningData.coarsestPartitionedHypergraph();

    if ( !is_vcycle ) {
      DegreeZeroHypernodeRemover degree_zero_hn_remover(context);
      if ( context.initial_partitioning.remove_degree_zero_hns_before_ip ) {
        degree_zero_hn_remover.removeDegreeZeroHypernodes(phg.hypergraph());
      }

      Context ip_context(context);
      ip_context.refinement = context.initial_partitioning.refinement;
      disableTimerAndStats(context);
      switch ( context.initial_partitioning.mode ) {
        case Mode::direct:
          // The pool initial partitioner consist of several flat bipartitioning
          // techniques. This case runs as a base case (k = 2) within recursive bipartitioning
          // or the deep multilevel scheme.
          pool::bipartition(phg, ip_context); break;
        // k-way partitions can be computed either by recursive bipartitioning
        // or deep multilevel partitioning.
        case Mode::recursive_bipartitioning:
          recursive_bipartitioning::partition(phg, ip_context); break;
        case Mode::deep_multilevel:
          deep_multilevel::partition(phg, ip_context); break;
        case Mode::UNDEFINED: ERROR("Undefined initial partitioning algorithm");
      }
      degree_zero_hn_remover.restoreDegreeZeroHypernodes(phg);
    } else {
      // When performing a V-cycle, we store the block IDs
      // of the input hypergraph as community IDs
      const Hypergraph& hypergraph = phg.hypergraph();
      phg.doParallelForAllNodes([&](const HypernodeID hn) {
        const PartitionID part_id = hypergraph.communityID(hn);
        ASSERT(part_id != kInvalidPartition && part_id < context.partition.k);
        ASSERT(phg.partID(hn) == kInvalidPartition);
        phg.setOnlyNodePart(hn, part_id);
      });
      phg.initializePartition();
    }

    io::printPartitioningResults(phg, context, "Initial Partitioning Results:");
    if ( context.partition.verbose_output && !is_vcycle ) {
      utils::Utilities::instance().getInitialPartitioningStats(
        context.utility_id).printInitialPartitioningStats();
    }
    timer.stop_timer("initial_partitioning");

    // ################## UNCOARSENING ##################
    io::printLocalSearchBanner(context);
    timer.start_timer("refinement", "Refinement");
    std::unique_ptr<IRefiner> label_propagation =
      LabelPropagationFactory::getInstance().createObject(
        context.refinement.label_propagation.algorithm, hypergraph, context);
    std::unique_ptr<IRefiner> fm =
      FMFactory::getInstance().createObject(
        context.refinement.fm.algorithm, hypergraph, context);

    std::unique_ptr<IUncoarsener> uncoarsener(nullptr);
    if (uncoarseningData.nlevel) {
      uncoarsener = std::make_unique<NLevelUncoarsener>(hypergraph, context, uncoarseningData);
    } else {
      uncoarsener = std::make_unique<MultilevelUncoarsener>(hypergraph, context, uncoarseningData);
    }
    partitioned_hg = uncoarsener->uncoarsen(label_propagation, fm);

    io::printPartitioningResults(partitioned_hg, context, "Local Search Results:");
    timer.stop_timer("refinement");

    return partitioned_hg;
  }
}

PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context) {
  PartitionedHypergraph partitioned_hg = multilevel_partitioning(hypergraph, context, false);

  // ################## V-CYCLES ##################
  if ( context.partition.num_vcycles > 0 && context.type == ContextType::main ) {
    partitionVCycle(hypergraph, partitioned_hg, context);
  }

  return partitioned_hg;
}

void partitionVCycle(Hypergraph& hypergraph,
                     PartitionedHypergraph& partitioned_hg,
                     const Context& context) {
  ASSERT(context.partition.num_vcycles > 0);

  for ( size_t i = 0; i < context.partition.num_vcycles; ++i ) {
    // Reset memory pool
    hypergraph.reset();
    parallel::MemoryPool::instance().reset();
    parallel::MemoryPool::instance().release_mem_group("Preprocessing");

    if ( context.partition.paradigm == Paradigm::nlevel ) {
      // Workaround: reset() function of hypergraph reinserts all removed hyperedges again.
      LargeHyperedgeRemover large_he_remover(context);
      large_he_remover.removeLargeHyperedgesInNLevelVCycle(hypergraph);
    }

    // The block IDs of the current partition are stored as community IDs.
    // This way coarsening does not contract nodes that do not belong to same block
    // of the input partition. For initial partitioning, we use the community IDs of
    // smallest hypergraph as initial partition.
    hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
      hypergraph.setCommunityID(hn, partitioned_hg.partID(hn));
    });

    // Perform V-cycle
    io::printVCycleBanner(context, i + 1);
    partitioned_hg = multilevel_partitioning(hypergraph, context, true /* V-cycle flag */ );
  }
}

}

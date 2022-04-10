/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/multilevel.h"

#include <memory>

#include "tbb/task.h"

#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/initial_partitioning/flat/pool_initial_partitioner.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/coarsening/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"
#include "mt-kahypar/partition/initial_partitioning/greedy_judicious_initial_partitioner.h"

namespace mt_kahypar::multilevel {

  class RefinementTask : public tbb::task {

  public:
    RefinementTask(Hypergraph& hypergraph,
                   PartitionedHypergraph& partitioned_hypergraph,
                   const Context& context,
                   std::shared_ptr<UncoarseningData> uncoarseningData) :
            _sparsifier(nullptr),
            _ip_context(context),
            _degree_zero_hn_remover(context),
            _uncoarsener(nullptr),
            _hg(hypergraph),
            _partitioned_hg(partitioned_hypergraph),
            _context(context),
            _uncoarseningData(uncoarseningData) {
      // Must be empty, because final partitioned hypergraph
      // is moved into this object
      _sparsifier = HypergraphSparsifierFactory::getInstance().createObject(
              _context.sparsification.similiar_net_combiner_strategy, _context);

      // Switch refinement context from IP to main
      _ip_context.refinement = _context.initial_partitioning.refinement;
      }

    tbb::task* execute() override {
      enableTimerAndStats();

      if ( _sparsifier->isSparsified() ) {
        // In that case, the sparsified hypergraph generated by the
        // heavy hyperedge remover was used for initial partitioning.
        // => Partition has to mapped from sparsified hypergraph to
        // coarsest partitioned hypergraph.
        io::printPartitioningResults(_sparsifier->sparsifiedPartitionedHypergraph(),
                                     _context, "Sparsified Initial Partitioning Results:");
        _degree_zero_hn_remover.restoreDegreeZeroHypernodes(
          _sparsifier->sparsifiedPartitionedHypergraph());
        _sparsifier->undoSparsification(_uncoarseningData->coarsestPartitionedHypergraph());
      } else {
        _degree_zero_hn_remover.restoreDegreeZeroHypernodes(
          _uncoarseningData->coarsestPartitionedHypergraph());
      }

      utils::Timer::instance().stop_timer("initial_partitioning");

      io::printPartitioningResults(_uncoarseningData->coarsestPartitionedHypergraph(),
                                   _context, "Initial Partitioning Results:");
      if ( _context.partition.verbose_output ) {
        utils::InitialPartitioningStats::instance().printInitialPartitioningStats();
      }

      // ################## LOCAL SEARCH ##################
      io::printLocalSearchBanner(_context);

      utils::Timer::instance().start_timer("refinement", "Refinement");
      std::unique_ptr<IRefiner> label_propagation =
              LabelPropagationFactory::getInstance().createObject(
                      _context.refinement.label_propagation.algorithm,
                      _hg, _context);
      std::unique_ptr<IRefiner> fm =
              FMFactory::getInstance().createObject(
                      _context.refinement.fm.algorithm,
                      _hg, _context);

      if (_uncoarseningData->nlevel) {
        _uncoarsener = std::make_unique<NLevelUncoarsener>(_hg, _context, *_uncoarseningData);
      } else {
        _uncoarsener = std::make_unique<MultilevelUncoarsener>(_hg, _context, *_uncoarseningData);
      }
      _partitioned_hg = _uncoarsener->uncoarsen(label_propagation, fm);
      utils::Timer::instance().stop_timer("refinement");

      io::printPartitioningResults(_partitioned_hg, _context, "Local Search Results:");

      return nullptr;
    }

  public:
    std::unique_ptr<IHypergraphSparsifier> _sparsifier;
    Context _ip_context;
    DegreeZeroHypernodeRemover _degree_zero_hn_remover;

  private:
    void enableTimerAndStats() {
      if ( _context.type == kahypar::ContextType::main ) {
        parallel::MemoryPool::instance().activate_unused_memory_allocations();
        utils::Timer::instance().enable();
        utils::Stats::instance().enable();
      }
    }

    std::unique_ptr<IUncoarsener> _uncoarsener;
    Hypergraph& _hg;
    PartitionedHypergraph& _partitioned_hg;
    const Context& _context;
    std::shared_ptr<UncoarseningData> _uncoarseningData;
  };

  class CoarseningTask : public tbb::task {

  public:
    CoarseningTask(Hypergraph& hypergraph,
                   IHypergraphSparsifier& sparsifier,
                   const Context& context,
                   const Context& ip_context,
                   DegreeZeroHypernodeRemover& degree_zero_hn_remover,
                   const bool vcycle,
                   UncoarseningData& uncoarseningData) :
            _hg(hypergraph),
            _sparsifier(sparsifier),
            _context(context),
            _ip_context(ip_context),
            _degree_zero_hn_remover(degree_zero_hn_remover),
            _vcycle(vcycle),
            _uncoarseningData(uncoarseningData) { }

    tbb::task* execute() override {
      // ################## COARSENING ##################

      if (_context.coarsening.skip_coarsening) {
        _uncoarseningData.finalizeCoarsening();
        utils::Timer::instance().start_timer("initial_partitioning", "Initial Partitioning");
        initialPartition(_uncoarseningData.coarsestPartitionedHypergraph());
        return nullptr;
      }
      mt_kahypar::io::printCoarseningBanner(_context);
      utils::Timer::instance().start_timer("coarsening", "Coarsening");
      _coarsener = CoarsenerFactory::getInstance().createObject(
              _context.coarsening.algorithm, _hg, _context, _uncoarseningData);
      _coarsener->coarsen();
      utils::Timer::instance().stop_timer("coarsening");

      Hypergraph& coarsestHypergraph = _coarsener->coarsestHypergraph();
      _coarsener.reset();
      if (_context.partition.verbose_output) {
        mt_kahypar::io::printHypergraphInfo(
                coarsestHypergraph, "Coarsened Hypergraph",
                _context.partition.show_memory_consumption);
      }

      // ################## INITIAL PARTITIONING ##################
      utils::Timer::instance().start_timer("initial_partitioning", "Initial Partitioning");
      if ( _context.useSparsification() ) {
        // Sparsify Hypergraph, if heavy hyperedge removal is enabled
        utils::Timer::instance().start_timer("sparsify_hypergraph", "Sparsify Hypergraph");
        _sparsifier.sparsify(coarsestHypergraph);
        utils::Timer::instance().stop_timer("sparsify_hypergraph");
      }

      if ( _sparsifier.isSparsified() ) {
        if (_context.partition.verbose_output) {
          mt_kahypar::io::printHypergraphInfo(
                  _sparsifier.sparsifiedHypergraph(), "Sparsified Hypergraph",
                  _context.partition.show_memory_consumption);
        }
        initialPartition(_sparsifier.sparsifiedPartitionedHypergraph());
      } else {
        initialPartition(_uncoarseningData.coarsestPartitionedHypergraph());
      }

      return nullptr;
    }

  private:
    void initialPartition(PartitionedHypergraph& phg) {
      io::printInitialPartitioningBanner(_context);

      if (_context.initial_partitioning.refinement.judicious.use_judicious_refinement) {
        judiciousIP(phg);
      } else if ( !_vcycle ) {
        if ( _context.initial_partitioning.remove_degree_zero_hns_before_ip ) {
          _degree_zero_hn_remover.removeDegreeZeroHypernodes(phg.hypergraph());
        }

        if ( _context.initial_partitioning.mode == InitialPartitioningMode::direct ) {
          disableTimerAndStats();
          PoolInitialPartitionerContinuation& ip_continuation = *new(allocate_continuation())
                  PoolInitialPartitionerContinuation(phg, _ip_context);
          spawn_initial_partitioner(ip_continuation);
        } else {
          std::unique_ptr<IInitialPartitioner> initial_partitioner =
                  InitialPartitionerFactory::getInstance().createObject(
                          _ip_context.initial_partitioning.mode, phg, _ip_context);
          initial_partitioner->initialPartition();
        }
      } else {
        // V-Cycle: Partition IDs are given by its community IDs
        const Hypergraph& hypergraph = phg.hypergraph();
        phg.doParallelForAllNodes([&](const HypernodeID hn) {
          const PartitionID part_id = hypergraph.communityID(hn);
          ASSERT(part_id != kInvalidPartition && part_id < _context.partition.k);
          ASSERT(phg.partID(hn) == kInvalidPartition);
          phg.setOnlyNodePart(hn, part_id);
        });
        phg.initializePartition();
      }
    }

    void disableTimerAndStats() {
      if ( _context.type == kahypar::ContextType::main ) {
        parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
        utils::Timer::instance().disable();
        utils::Stats::instance().disable();
      }
    }

    void judiciousIP(PartitionedHypergraph& phg) {
      GreedyJudiciousInitialPartitioner ip(phg, _ip_context);
      ip.initialPartition(_context.partition.seed);
      const size_t num_runs = _context.initial_partitioning.runs;
      std::uniform_int_distribution<> distrib(0, std::numeric_limits<int>::max());
      std::mt19937 g(_context.partition.seed);
      tbb::task_group tg;
      vec<std::pair<HyperedgeWeight, vec<size_t>>> partitions(num_runs);
      tbb::enumerable_thread_specific<PartitionedHypergraph> phgs([&] { return construct_phg(phg); });
      auto ip_run = [&](const size_t seed, const size_t i) {
        auto &local_phg = phgs.local();
        // run IP and extract part IDs
        GreedyJudiciousInitialPartitioner ip(local_phg, _ip_context);
        ip.initialPartition(seed);
        partitions[i].second.resize(phg.initialNumNodes());
        for (size_t j = 0; j < phg.initialNumNodes(); ++j) {
          partitions[i].second[j] = local_phg.partID(j);
        }
        partitions[i].first = metrics::judiciousLoad(local_phg);
      };
      for(size_t i = 0; i < num_runs; ++i) {
        tg.run(std::bind(ip_run, distrib(g), i));
      }
      tg.wait();
      phg.resetPartition();
      // choose best partititon and assign it to the hypergraph
      auto& best_partition = *std::min_element(partitions.begin(), partitions.end());
      for (size_t i = 0; i < phg.initialNumNodes(); ++i) {
        phg.setNodePart(i, best_partition.second[i]);
      }
    }

    PartitionedHypergraph construct_phg(PartitionedHypergraph& phg) {
      return PartitionedHypergraph(_context.partition.k, phg.hypergraph());
    }

    Hypergraph& _hg;
    IHypergraphSparsifier& _sparsifier;
    const Context& _context;
    const Context& _ip_context;
    DegreeZeroHypernodeRemover& _degree_zero_hn_remover;
    std::unique_ptr<ICoarsener> _coarsener;
    const bool _vcycle;
    UncoarseningData& _uncoarseningData;
  };

// ! Helper function that spawns the multilevel partitioner in
// ! TBB continuation style with a given parent task.
  static void spawn_multilevel_partitioner(Hypergraph& hypergraph,
                                           PartitionedHypergraph& partitioned_hypergraph,
                                           const Context& context,
                                           const bool vcycle,
                                           tbb::task& parent) {
    // The coarsening task is first executed and once it finishes the
    // refinement task continues (without blocking)
    bool nlevel = context.coarsening.algorithm == CoarseningAlgorithm::nlevel_coarsener;
    std::shared_ptr<UncoarseningData> uncoarseningData =
      std::make_shared<UncoarseningData>(nlevel, hypergraph, context);

    RefinementTask& refinement_task = *new(parent.allocate_continuation())
            RefinementTask(hypergraph, partitioned_hypergraph, context, uncoarseningData);
    refinement_task.set_ref_count(1);
    CoarseningTask& coarsening_task = *new(refinement_task.allocate_child()) CoarseningTask(
            hypergraph, *refinement_task._sparsifier, context, refinement_task._ip_context,
            refinement_task._degree_zero_hn_remover, vcycle, *uncoarseningData);
    tbb::task::spawn(coarsening_task);
  }

  class MultilevelPartitioningTask : public tbb::task {

  public:
    MultilevelPartitioningTask(Hypergraph& hypergraph,
                               PartitionedHypergraph& partitioned_hypergraph,
                               const Context& context,
                               const bool vcycle) :
            _hg(hypergraph),
            _partitioned_hg(partitioned_hypergraph),
            _context(context),
            _vcycle(vcycle) { }

    tbb::task* execute() override {
      spawn_multilevel_partitioner(
              _hg, _partitioned_hg, _context, _vcycle, *this);
      return nullptr;
    }

  private:
    Hypergraph& _hg;
    PartitionedHypergraph& _partitioned_hg;
    const Context& _context;
    const bool _vcycle;
  };


PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context, const bool vcycle) {
  PartitionedHypergraph partitioned_hypergraph;
  MultilevelPartitioningTask& multilevel_task = *new(tbb::task::allocate_root())
          MultilevelPartitioningTask(hypergraph, partitioned_hypergraph, context, vcycle);
  tbb::task::spawn_root_and_wait(multilevel_task);
  return partitioned_hypergraph;
}



void partition_async(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context, tbb::task* parent) {
  ASSERT(parent);
  spawn_multilevel_partitioner(hypergraph, partitioned_hypergraph, context, false, *parent);
}

}

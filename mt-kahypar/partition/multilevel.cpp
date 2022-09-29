/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/multilevel.h"

#include <memory>

#include "tbb/task.h"

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/partition/coarsening/separated_nodes/snodes_sync_coarsening.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/initial_partitioning/flat/pool_initial_partitioner.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/coarsening/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"

namespace mt_kahypar::multilevel {
  using ds::Array;

  void reconstructHierarchyWithSeparatedNodes(Hypergraph& hg, const Context& context,
                                              SepNodesStack& stack, UncoarseningData& uncoarseningData,
                                              bool with_partition=true) {
      utils::Timer::instance().start_timer("reconstruct_hierarchy", "Reconstruct Hierarchy");

      vec<Level>& levels = uncoarseningData.hierarchy;
      HypernodeWeight correct_weight = hg.totalActiveWeight();
      unused(correct_weight);
      ASSERT([&] {
        correct_weight += hg.separatedNodes().finest().initialWeight();
        return true;
      }());
      ASSERT(stack.numLevels() == levels.size() + 1);

      for (size_t i = 0; i < levels.size(); ++i) {
        Hypergraph& old_hg = (i == 0) ? hg : levels[i - 1].contractedHypergraph();
        const HypernodeID num_graph_nodes = levels[i].contractedHypergraph().initialNumNodes();
        const vec<HypernodeID>& old_communities = levels[i].communities();
        const vec<HypernodeID>& sep_mapping = stack.mapping(i, false);
        SeparatedNodes& s_nodes = stack.atLevel(i, false);
        ASSERT(s_nodes.numGraphNodes() == old_hg.initialNumNodes());
        ASSERT(old_communities.size() == old_hg.initialNumNodes());
        ASSERT(sep_mapping.size() >= s_nodes.numNodes());

        const size_t old_size = old_communities.size();
        vec<HypernodeID> new_communities(old_size + s_nodes.currentBatchIndex());
        tbb::parallel_for(0UL, new_communities.size(), [&](const size_t& pos) {
          if (pos < old_size) {
            ASSERT(old_communities[pos] == kInvalidHypernode || old_communities[pos] < num_graph_nodes);
            new_communities[pos] = old_communities[pos];
          } else {
            ASSERT(sep_mapping[pos - old_size] < stack.atLevel(i + 1, false).numNodes());
            new_communities[pos] = sep_mapping[pos - old_size] + num_graph_nodes;
          }
        });
        levels[i].communities() = std::move(new_communities);

        vec<HypernodeID>& communities = levels[i].communities();
        tbb::parallel_for(s_nodes.currentBatchIndex(), s_nodes.numNodes(), [&](const HypernodeID& node) {
          const HypernodeID original_id = s_nodes.originalHypernodeID(node);
          ASSERT(original_id != kInvalidHypernode && communities[original_id] == kInvalidHypernode);
          communities[original_id] = sep_mapping[node] + num_graph_nodes;
        });

        s_nodes.popBatch();
        if (s_nodes.numNodes() > 0) {
          old_hg = HypergraphFactory::reinsertSeparatedNodes(old_hg, s_nodes);
        }
        old_hg.setSeparatedNodes(nullptr);
        ASSERT(correct_weight == old_hg.totalWeight(), V(correct_weight) << V(old_hg.totalWeight()));
      }

      Hypergraph& last_hg = (levels.size() == 0) ? hg : levels.back().contractedHypergraph();
      PartitionedHypergraph& phg = *uncoarseningData.partitioned_hg;
      Array<PartIdType> part_ids;
      if (with_partition) {
        part_ids.resize(phg.partIDsSize(context.type != kahypar::ContextType::main), PartIdType(kInvalidPartition));
        phg.extractPartIDs(part_ids, false);
      }

      if (stack.coarsest().numNodes() > 0) {
        last_hg = HypergraphFactory::reinsertSeparatedNodes(last_hg, stack.coarsest());
      }
      last_hg.setSeparatedNodes(nullptr);
      ASSERT(correct_weight == last_hg.totalWeight());
      phg.setHypergraph(last_hg);

      if (with_partition) {
        phg.resetData();
        phg.doParallelForAllNodes([&](const HypernodeID node) {
          phg.setOnlyNodePart(node, part_ids[node].load());
        });
        phg.initializePartition();
      }

      ASSERT([&] {
        for (size_t i = 0; i < levels.size(); ++i) {
          const Hypergraph& old_hg = (i == 0) ? hg : levels[i - 1].contractedHypergraph();
          if (old_hg.initialNumNodes() != levels[i].communities().size()) {
            return false;
          }
          for (const HypernodeID& node: levels[i].communities()) {
            if (node == kInvalidHypernode || node >= levels[i].contractedHypergraph().initialNumNodes()) {
              return false;
            }
          }
        }
        return true;
      }());
      ASSERT(!with_partition || [&] {
        for (const HypernodeID& node: phg.nodes()) {
          if (phg.partID(node) == kInvalidPartition) {
            return false;
          }
        }
        return true;
      }());

      utils::Timer::instance().stop_timer("reconstruct_hierarchy");
    }

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

      // undo replacement of IP hypergraph for separated nodes
      _uncoarseningData->initializeRefinement();

      utils::Timer::instance().stop_timer("initial_partitioning");

      io::printPartitioningResults(_uncoarseningData->coarsestPartitionedHypergraph(),
                                   _context, "Initial Partitioning Results:");
      if (_context.graphviz_file != "" && _context.type == kahypar::ContextType::main) {
        #ifdef USE_GRAPH_PARTITIONER
        utils::outputGraphvizFile(_uncoarseningData->coarsestPartitionedHypergraph(), _context.graphviz_file, false, ".result_plain");
        utils::outputGraphvizFile(_uncoarseningData->coarsestPartitionedHypergraph(), _context.graphviz_file, true, ".result_incl");
        #endif
      }
      if ( _context.partition.verbose_output ) {
        utils::InitialPartitioningStats::instance().printInitialPartitioningStats();
      }

      // ################## SEPARATED NODES COARSENING ##################
      if (_context.refinement.include_separated) {
        Array<PartitionID> part_ids;
        Array<PartitionID>* part_id_ptr = nullptr;
        if (_context.refinement.separated_partition_aware_coarsening) {
          _uncoarseningData->partitioned_hg->resetSeparatedParts();
          const SeparatedNodes& s_nodes = _hg.separatedNodes().finest();
          const HyperedgeWeight added_cut = star_partitioning::partition(*_uncoarseningData->partitioned_hg,
                                                                         s_nodes, _context, false);
          part_ids.resize(s_nodes.numNodes(), kInvalidPartition);
          tbb::parallel_for(ID(0), s_nodes.numNodes(), [&] (const HypernodeID& node) {
            part_ids[node] = _uncoarseningData->partitioned_hg->separatedPartID(node);
          });
          part_id_ptr = &part_ids;
        }
        SepNodesStack stack(_hg.separatedNodes().finest().createCopyFromSavepoint());
        const HypernodeID start_num_nodes = _hg.numSeparatedNodes();
        const HypernodeID target_num_nodes = _uncoarseningData->calculateSeparatedNodesTargetSize(
                                              _uncoarseningData->coarsestPartitionedHypergraph().hypergraph(), _hg);
        star_partitioning::coarsenSynchronized(stack, _hg, _uncoarseningData->hierarchy, _context,
                                               start_num_nodes, target_num_nodes, part_id_ptr);
        _uncoarseningData->partitioned_hg->setSeparatedNodes(&stack);
        if (_context.refinement.separated_partition_aware_coarsening) {
          vec<CAtomic<PartitionID>>& graph_part_ids = _uncoarseningData->partitioned_hg->separatedPartIDs();
          ASSERT(graph_part_ids.size() >= part_ids.size());
          tbb::parallel_for(0UL, graph_part_ids.size(), [&] (const size_t& pos) {
            if (pos < part_ids.size()) {
              graph_part_ids[pos].store(part_ids[pos]);
            } else {
              graph_part_ids[pos].store(kInvalidPartition);
            }
          });
        } else {
          // TODO!!
          _uncoarseningData->partitioned_hg->resetSeparatedParts();
          const HyperedgeWeight added_cut = star_partitioning::partition(*_uncoarseningData->partitioned_hg,
                                                                         stack.coarsest(), _context, false);
        }

        reconstructHierarchyWithSeparatedNodes(_hg, _context, stack, *_uncoarseningData);
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
      if ( _context.type == kahypar::ContextType::main && _context.partition.mode == Mode::direct ) {
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
      if (_context.refinement.include_separated || _context.coarsening.sep_nodes_sync_coarsening) {
        // separated nodes
        _hg.separatedNodes().finest().setSavepoint();
      }

      // ################## COARSENING ##################
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
      if (_context.coarsened_stats_file != "") {
        utils::printDistributionStatsToCSV(coarsestHypergraph, _context.coarsened_stats_file);
      }
      if (_context.graphviz_file != "") {
        #ifdef USE_GRAPH_PARTITIONER
        utils::outputGraphvizFile(coarsestHypergraph, _context.graphviz_file, false, ".coarse_plain");
        utils::outputGraphvizFile(coarsestHypergraph, _context.graphviz_file, true, ".coarse_incl");
        #endif
      }

      // ################## SEPARATED NODES COARSENING ##################
      if (_context.coarsening.sep_nodes_sync_coarsening) {
        SepNodesStack stack(_hg.separatedNodes().finest().createCopyFromSavepoint());
        const HypernodeID start_num_nodes = _hg.numSeparatedNodes();
        const HypernodeID target_num_nodes = _uncoarseningData.calculateSeparatedNodesTargetSize(
                                              _uncoarseningData.coarsestPartitionedHypergraph().hypergraph(), _hg);
        star_partitioning::coarsenSynchronized(stack, _hg, _uncoarseningData.hierarchy, _context,
                                               start_num_nodes, target_num_nodes, nullptr);
        _uncoarseningData.partitioned_hg->setSeparatedNodes(&stack);
        _uncoarseningData.partitioned_hg->resetSeparatedParts();
        reconstructHierarchyWithSeparatedNodes(_hg, _context, stack, _uncoarseningData, false);
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

      if ( !_vcycle ) {
        if ( _context.initial_partitioning.remove_degree_zero_hns_before_ip ) {
          _degree_zero_hn_remover.removeDegreeZeroHypernodes(phg.hypergraph());
        }

        if ( _context.initial_partitioning.mode == Mode::direct ) {
          if (_context.initial_partitioning.rater == IPSNodesRater::tracker) {
            phg.separatedNodes().finest().initializeOutwardEdges();
          }
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
      if ( _context.type == kahypar::ContextType::main && _context.partition.mode == Mode::direct ) {
        parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
        utils::Timer::instance().disable();
        utils::Stats::instance().disable();
      }
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

  class VCycleTask : public tbb::task {
  public:
    VCycleTask(Hypergraph& hypergraph,
               PartitionedHypergraph& partitioned_hypergraph,
               const Context& context) :
            _hg(hypergraph),
            _partitioned_hg(partitioned_hypergraph),
            _context(context) { }

    tbb::task* execute() override {
      ASSERT(_context.partition.num_vcycles > 0);

      for ( size_t i = 0; i < _context.partition.num_vcycles; ++i ) {
        // Reset memory pool
        _hg.reset();
        parallel::MemoryPool::instance().reset();
        parallel::MemoryPool::instance().release_mem_group("Preprocessing");

        if ( _context.partition.paradigm == Paradigm::nlevel ) {
          // Workaround: reset() function of hypergraph reinserts all removed
          // hyperedges to incident net lists of each vertex again.
          LargeHyperedgeRemover large_he_remover(_context);
          large_he_remover.removeLargeHyperedgesInNLevelVCycle(_hg);
        }

        // Store partition and assign it as community ids in order to
        // restrict contractions in v-cycle to partition ids
        _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          _hg.setCommunityID(hn, _partitioned_hg.partID(hn));
        });

        // V-Cycle Multilevel Partitioning
        io::printVCycleBanner(_context, i + 1);
        MultilevelPartitioningTask& multilevel_task = *new(tbb::task::allocate_root())
                MultilevelPartitioningTask(_hg, _partitioned_hg, _context, true /* vcycle */);
        tbb::task::spawn_root_and_wait(multilevel_task);
      }

      return nullptr;
    }

  private:
    Hypergraph& _hg;
    PartitionedHypergraph& _partitioned_hg;
    const Context& _context;
  };


PartitionedHypergraph partition(Hypergraph& hypergraph, const Context& context) {
  PartitionedHypergraph partitioned_hypergraph;
    MultilevelPartitioningTask& multilevel_task = *new(tbb::task::allocate_root())
            MultilevelPartitioningTask(hypergraph, partitioned_hypergraph, context, false);
    tbb::task::spawn_root_and_wait(multilevel_task);

    if ( context.partition.num_vcycles > 0 && context.type == kahypar::ContextType::main ) {
      partitionVCycle(hypergraph, partitioned_hypergraph, context);
    }
  return partitioned_hypergraph;
}


void partition_async(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context, tbb::task* parent) {
  ASSERT(parent);

  if ( context.partition.num_vcycles > 0 && context.type == kahypar::ContextType::main ) {
    VCycleTask& vcycle_task = *new(parent->allocate_continuation())
            VCycleTask(hypergraph, partitioned_hypergraph, context);
    MultilevelPartitioningTask& multilevel_task = *new(vcycle_task.allocate_child())
            MultilevelPartitioningTask(hypergraph, partitioned_hypergraph, context, false);
    tbb::task::spawn(multilevel_task);
  } else {
    spawn_multilevel_partitioner(hypergraph, partitioned_hypergraph, context, false, *parent);
  }
}


void partitionVCycle(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                     const Context& context) {
  VCycleTask& vcycle_task = *new(tbb::task::allocate_root())
          VCycleTask(hypergraph, partitioned_hypergraph, context);
  tbb::task::spawn_root_and_wait(vcycle_task);
}

}

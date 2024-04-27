/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include <mt-kahypar/partition/coarsening/multilevel_uncoarsener.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/cast.h"
#include <fstream>

namespace mt_kahypar {

  template<typename TypeTraits>
  void MultilevelUncoarsener<TypeTraits>::initializeImpl() {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    _current_metrics = Base::initializeMetrics(partitioned_hg);
    Base::initializeRefinementAlgorithms();

    if (_context.type == ContextType::main) {
      _context.initial_km1 = _current_metrics.quality;
    }

    // Enable progress bar if verbose output is enabled
    if ( _context.partition.verbose_output && _context.partition.enable_progress_bar && !debug ) {
      _progress.enable();
      _progress.setObjective(_current_metrics.quality);
    }

    // Pass target graph to partitioned hypergraph
    partitioned_hg.setTargetGraph(_target_graph);

    _current_level = _uncoarseningData.hierarchy.size();
    _num_levels = _current_level;
  }

  template<typename TypeTraits>
  bool MultilevelUncoarsener<TypeTraits>::isTopLevelImpl() const {
    return _current_level < 0;
  }

  template<typename TypeTraits>
  void MultilevelUncoarsener<TypeTraits>::projectToNextLevelAndRefineImpl() {
    PartitionedHypergraph& partitioned_hg = *_uncoarseningData.partitioned_hg;
    std::cout << "test\n\n";
    if ( _current_level == _num_levels) {
      // We always start with a refinement pass on the smallest hypergraph.
      // The next calls to this function will then project the partition to the next level
      // and perform refinement until we reach the input hypergraph.
      IUncoarsener<TypeTraits>::refine();
      _progress.setObjective(_current_metrics.quality);
      _progress += partitioned_hg.initialNumNodes();
    } else {
      ASSERT(_current_level >= 0);
      // Project partition to the hypergraph on the next level of the hierarchy
      _timer.start_timer("projecting_partition", "Projecting Partition");
      const size_t num_nodes_on_previous_level = partitioned_hg.initialNumNodes();
      if (_current_level == 0) {
        partitioned_hg.setHypergraph(_hg);
      } else {
        partitioned_hg.setHypergraph((_uncoarseningData.hierarchy)[_current_level-1].contractedHypergraph());
      }
      // Hypergraph stores partition from previous level.
      // We now extract the block ids and reset the partition to reuse the
      // data structure for the current level.
      partitioned_hg.extractPartIDs(_block_ids);
      partitioned_hg.resetData();
      GainCachePtr::resetGainCache(_gain_cache);


      /*readPartitionFile(partitioned_hg, _context.partition.partition_input_file);*/

      // Assign nodes of current level to their corresponding representative of the previous level
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const HypernodeID coarse_hn = (_uncoarseningData.hierarchy)[_current_level].mapToContractedHypergraph(hn);
        const PartitionID block = _block_ids[coarse_hn];
        ASSERT(block != kInvalidPartition && block < partitioned_hg.k());
        partitioned_hg.setOnlyNodePart(hn, block);
      });
      if(_current_level==0){
        std::ifstream myfile; 
        myfile.open(_context.partition.partition_input_file);
        HypernodeID hn = 0;
        while(myfile){
          std::string nextline;
          std::getline(myfile, nextline);
          std::cout << nextline;
          if(nextline != ""){
            ASSERT(partitioned_hg.k() == 2);
            ASSERT(stoi(nextline) < partitioned_hg.k());
            partitioned_hg.setOnlyNodePart(hn, stoi(nextline));
          }        
        }
        myfile.close();
      }
      partitioned_hg.initializePartition();
      _timer.stop_timer("projecting_partition");

      std::cout << "quality: " << metrics::quality(partitioned_hg, _context) << "\n\n\n";

      // Improve partition
      IUncoarsener<TypeTraits>::refine();

      // Update Progress Bar
      _progress.setObjective(_current_metrics.quality);
      _progress += partitioned_hg.initialNumNodes() - num_nodes_on_previous_level;
      std::cout << "done\n";
    }

    /*ASSERT(metrics::quality(*_uncoarseningData.partitioned_hg, _context) == _current_metrics.quality,
      V(_current_metrics.quality) << V(metrics::quality(*_uncoarseningData.partitioned_hg, _context)));*/

    --_current_level;
  }

  template<typename TypeTraits>
  void MultilevelUncoarsener<TypeTraits>::rebalancingImpl() {
    // If we reach the top-level hypergraph and the partition is still imbalanced,
    // we use a rebalancing algorithm to restore balance.
    std::cout << static_cast<int>(_context.type == ContextType::main);
    if (_context.type == ContextType::main && !metrics::isBalanced(*_uncoarseningData.partitioned_hg, _context)) {
      std::cout << "rebalancing\n\n\n";
      const HyperedgeWeight quality_before = _current_metrics.quality;
      if (_context.partition.verbose_output) {
        LOG << RED << "Partition is imbalanced (Current Imbalance:"
        << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;

        LOG << "Part weights: (violations in red)";
        io::printPartWeightsAndSizes(*_uncoarseningData.partitioned_hg, _context);
      }

      if ( !_context.partition.deterministic ) {
        if (_context.partition.verbose_output) {
          LOG << RED << "Start rebalancing!" << END;
        }

        // Preform rebalancing
        _timer.start_timer("rebalance", "Rebalance");
        mt_kahypar_partitioned_hypergraph_t phg =
          utils::partitioned_hg_cast(*_uncoarseningData.partitioned_hg);
        _rebalancer->refine(phg, {}, _current_metrics, 0.0);
        _timer.stop_timer("rebalance");

        const HyperedgeWeight quality_after = _current_metrics.quality;
        if (_context.partition.verbose_output) {
          const HyperedgeWeight quality_delta = quality_after - quality_before;
          if (quality_delta > 0) {
            LOG << RED << "Rebalancer decreased solution quality by" << quality_delta
            << "(Current Imbalance:" << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;
          } else {
            LOG << GREEN << "Rebalancer improves solution quality by" << abs(quality_delta)
            << "(Current Imbalance:" << metrics::imbalance(*_uncoarseningData.partitioned_hg, _context) << ")" << END;
          }
        }
      } else {
        if (_context.partition.verbose_output) {
          LOG << RED << "Skip rebalancing since deterministic mode is activated" << END;
        }
      }


      ASSERT(metrics::quality(*_uncoarseningData.partitioned_hg, _context) == _current_metrics.quality,
        V(_current_metrics.quality) << V(metrics::quality(*_uncoarseningData.partitioned_hg, _context)));
    }
  }

  template<typename TypeTraits>
  HyperedgeWeight MultilevelUncoarsener<TypeTraits>::getObjectiveImpl() const {
    return _current_metrics.quality;
  }

  template<typename TypeTraits>
  void MultilevelUncoarsener<TypeTraits>::updateMetricsImpl() {
    _current_metrics = Base::initializeMetrics(*_uncoarseningData.partitioned_hg);
    _progress.setObjective(_current_metrics.quality);
  }

  template<typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph& MultilevelUncoarsener<TypeTraits>::currentPartitionedHypergraphImpl() {
    return *_uncoarseningData.partitioned_hg;
  }

  template<typename TypeTraits>
  HypernodeID MultilevelUncoarsener<TypeTraits>::currentNumberOfNodesImpl() const {
    return _uncoarseningData.partitioned_hg->initialNumNodes();
  }

  template<typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph&& MultilevelUncoarsener<TypeTraits>::movePartitionedHypergraphImpl() {
    ASSERT(isTopLevelImpl());
    std::cout << "before\n";
    return std::move(*_uncoarseningData.partitioned_hg);
    std::cout << "after\n";
  }

  template<typename TypeTraits>
  void MultilevelUncoarsener<TypeTraits>::refineImpl() {
    PartitionedHypergraph& partitioned_hypergraph = *_uncoarseningData.partitioned_hg;
    const double time_limit = Base::refinementTimeLimit(_context, (_uncoarseningData.hierarchy)[_current_level].coarseningTime());
    std::cout << "test\n";
    

    if ( debug && _context.type == ContextType::main ) {
      io::printHypergraphInfo(partitioned_hypergraph.hypergraph(),
        _context, "Refinement Hypergraph", false);
      DBG << "Start Refinement -" << _context.partition.objective << "=" << _current_metrics.quality
      << ", imbalance = " << _current_metrics.imbalance;
    }

    parallel::scalable_vector<HypernodeID> dummy;
    bool improvement_found = true;
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hypergraph);
    while( improvement_found ) {
      improvement_found = false;
      const HyperedgeWeight metric_before = _current_metrics.quality;

      if ( _rebalancer && _context.refinement.rebalancer != RebalancingAlgorithm::do_nothing ) {
        _rebalancer->initialize(phg);
      }

      if ( _label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        _label_propagation->initialize(phg);
        _timer.stop_timer("initialize_lp_refiner");

        _timer.start_timer("label_propagation", "Label Propagation");
        improvement_found |= _label_propagation->refine(phg, dummy, _current_metrics, time_limit);
        _timer.stop_timer("label_propagation");
      }

      if ( _fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        _fm->initialize(phg);
        _timer.stop_timer("initialize_fm_refiner");

        _timer.start_timer("fm", "FM");
        improvement_found |= _fm->refine(phg, dummy, _current_metrics, time_limit);
        _timer.stop_timer("fm");
      }

      if ( _flows && _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
        _timer.start_timer("initialize_flow_scheduler", "Initialize Flow Scheduler");
        _flows->initialize(phg);
        _timer.stop_timer("initialize_flow_scheduler");

        _timer.start_timer("flow_refinement_scheduler", "Flow Refinement Scheduler");
        improvement_found |= _flows->refine(phg, dummy, _current_metrics, time_limit);
        _timer.stop_timer("flow_refinement_scheduler");
      }

      if ( _context.type == ContextType::main ) {
        /*ASSERT(_current_metrics.quality == metrics::quality(partitioned_hypergraph, _context),
          "Actual metric" << V(metrics::quality(partitioned_hypergraph, _context)) <<
          "does not match the metric updated by the refiners" << V(_current_metrics.quality));*/
      }

      const HyperedgeWeight metric_after = _current_metrics.quality;
      const double relative_improvement = 1.0 -
        static_cast<double>(metric_after) / metric_before;
      if ( !_context.refinement.refine_until_no_improvement ||
           relative_improvement <= _context.refinement.relative_improvement_threshold ) {
        break;
      }
    }

    if ( _context.type == ContextType::main) {
      DBG << "--------------------------------------------------\n";
    }
  }


  template<typename TypeTraits>
  void MultilevelUncoarsener<TypeTraits>::doLastRefineImpl(PartitionedHypergraph *partitioned_hypergraph) {
    const double time_limit = 1000000;
    /*partitioned_hypergraph.resetData();*/
    if ( debug && _context.type == ContextType::main ) {
      io::printHypergraphInfo(partitioned_hypergraph->hypergraph(),
        _context, "Refinement Hypergraph", false);
      DBG << "Start Refinement -" << _context.partition.objective << "=" << _current_metrics.quality
      << ", imbalance = " << _current_metrics.imbalance;
    }
    initializeImpl();
    parallel::scalable_vector<HypernodeID> dummy;
    bool improvement_found = true;
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(*partitioned_hypergraph);
    while( improvement_found ) {
      improvement_found = false;
      const HyperedgeWeight metric_before = _current_metrics.quality;

      if ( _rebalancer && _context.refinement.rebalancer != RebalancingAlgorithm::do_nothing ) {
        _rebalancer->initialize(phg);
      }

      /*if ( _label_propagation && _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        std::cout << "lp\n";
        _timer.start_timer("initialize_lp_refiner", "Initialize LP Refiner");
        _label_propagation->initialize(phg);
        _timer.stop_timer("initialize_lp_refiner");

        _timer.start_timer("label_propagation", "Label Propagation");
        improvement_found |= _label_propagation->refine(phg, dummy, _current_metrics, time_limit);
        _timer.stop_timer("label_propagation");
      }*/

      if ( _fm && _context.refinement.fm.algorithm != FMAlgorithm::do_nothing ) {
        //std::cout << "fm\n";
        _timer.start_timer("initialize_fm_refiner", "Initialize FM Refiner");
        _fm->initialize(phg);
        _timer.stop_timer("initialize_fm_refiner");

        _timer.start_timer("fm", "FM");
        improvement_found |= _fm->refine(phg, dummy, _current_metrics, time_limit);
        _timer.stop_timer("fm");
      }

      if ( _flows && _context.refinement.flows.algorithm != FlowAlgorithm::do_nothing ) {
        std::cout << "flows\n";
        _timer.start_timer("initialize_flow_scheduler", "Initialize Flow Scheduler");
        _flows->initialize(phg);
        _timer.stop_timer("initialize_flow_scheduler");

        _timer.start_timer("flow_refinement_scheduler", "Flow Refinement Scheduler");
        improvement_found |= _flows->refine(phg, dummy, _current_metrics, time_limit);
        _timer.stop_timer("flow_refinement_scheduler");
      }

      if ( _context.type == ContextType::main ) {
        /*ASSERT(_current_metrics.quality == metrics::quality(partitioned_hypergraph, _context),
          "Actual metric" << V(metrics::quality(partitioned_hypergraph, _context)) <<
          "does not match the metric updated by the refiners" << V(_current_metrics.quality));*/
      }

      const HyperedgeWeight metric_after = _current_metrics.quality;
      const double relative_improvement = 1.0 -
        static_cast<double>(metric_after) / metric_before;
      if ( !_context.refinement.refine_until_no_improvement ||
           relative_improvement <= _context.refinement.relative_improvement_threshold ) {
        break;
      }
    }

    if ( _context.type == ContextType::main) {
      DBG << "--------------------------------------------------\n";
    }
    
  }

  INSTANTIATE_CLASS_WITH_TYPE_TRAITS(MultilevelUncoarsener)

}

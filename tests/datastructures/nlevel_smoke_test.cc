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

#include "gmock/gmock.h"

#include <atomic>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph.h"
#include "mt-kahypar/datastructures/dynamic_hypergraph_factory.h"
#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"

#include "smoke_test_common.h"

namespace mt_kahypar {
namespace ds {

DynamicHypergraph simulateNLevel(DynamicHypergraph& hypergraph,
                                 DynamicPartitionedHypergraph& partitioned_hypergraph,
                                 const BatchVector& contraction_batches,
                                 const size_t batch_size,
                                 const bool parallel) {

  auto timer_key = [&](const std::string& key) {
    if ( parallel ) {
      return key + "_parallel";
    } else {
      return key;
    }
  };

  parallel::scalable_vector<parallel::scalable_vector<ParallelHyperedge>> removed_hyperedges;
  for ( size_t i = 0; i < contraction_batches.size(); ++i ) {
    utils::Timer::instance().start_timer(timer_key("contractions"), "Contractions");
    const parallel::scalable_vector<Memento>& contractions = contraction_batches[i];
    if ( parallel ) {
      tbb::parallel_for(0UL, contractions.size(), [&](const size_t j) {
        const Memento& memento = contractions[j];
        hypergraph.registerContraction(memento.u, memento.v);
        hypergraph.contract(memento.v);
      });
    } else {
      for ( size_t j = 0; j < contractions.size(); ++j ) {
        const Memento& memento = contractions[j];
        hypergraph.registerContraction(memento.u, memento.v);
        hypergraph.contract(memento.v);
      }
    }
    utils::Timer::instance().stop_timer(timer_key("contractions"));

    utils::Timer::instance().start_timer(timer_key("remove_parallel_nets"), "Parallel Net Detection");
    removed_hyperedges.emplace_back(hypergraph.removeSinglePinAndParallelHyperedges());
    utils::Timer::instance().stop_timer(timer_key("remove_parallel_nets"));
  }

  utils::Timer::instance().start_timer(timer_key("copy_coarsest_hypergraph"), "Copy Coarsest Hypergraph");
  DynamicHypergraph coarsest_hypergraph;
  if ( parallel ) {
    coarsest_hypergraph = hypergraph.copy(parallel_tag_t());
  } else {
    coarsest_hypergraph = hypergraph.copy();
  }
  utils::Timer::instance().stop_timer(timer_key("copy_coarsest_hypergraph"));


  utils::Timer::instance().start_timer(timer_key("initial_partition"), "Initial Partition");

  {
    utils::Timer::instance().start_timer(timer_key("compactify_hypergraph"), "Compactify Hypergraph");
    auto res = DynamicHypergraphFactory::compactify(hypergraph);
    DynamicHypergraph& compactified_hg = res.first;
    auto& hn_mapping = res.second;
    DynamicPartitionedHypergraph compactified_phg(
      partitioned_hypergraph.k(), compactified_hg, parallel_tag_t());
    utils::Timer::instance().stop_timer(timer_key("compactify_hypergraph"));

    utils::Timer::instance().start_timer(timer_key("generate_random_partition"), "Generate Random Partition");
    generateRandomPartition(compactified_phg);
    utils::Timer::instance().stop_timer(timer_key("generate_random_partition"));

    utils::Timer::instance().start_timer(timer_key("project_partition"), "Project Partition");
    partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      partitioned_hypergraph.setOnlyNodePart(hn, compactified_phg.partID(hn_mapping[hn]));
    });
    utils::Timer::instance().stop_timer(timer_key("project_partition"));
  }

  utils::Timer::instance().start_timer(timer_key("initialize_partition"), "Initialize Partition");
  partitioned_hypergraph.initializePartition();
  utils::Timer::instance().stop_timer(timer_key("initialize_partition"));

  utils::Timer::instance().start_timer(timer_key("initialize_gain_cache"), "Initialize Initialize Gain Cache");
  partitioned_hypergraph.initializeGainCache();
  utils::Timer::instance().stop_timer(timer_key("initialize_gain_cache"));

  utils::Timer::instance().stop_timer(timer_key("initial_partition"));

  utils::Timer::instance().start_timer(timer_key("create_batch_uncontraction_hierarchy"), "Create n-Level Hierarchy");
  const size_t tmp_batch_size = parallel ? batch_size : 1;
  auto versioned_batches = hypergraph.createBatchUncontractionHierarchy(tmp_batch_size);
  utils::Timer::instance().stop_timer(timer_key("create_batch_uncontraction_hierarchy"));

  utils::Timer::instance().start_timer(timer_key("batch_uncontractions"), "Batch Uncontractions");
  while ( !versioned_batches.empty() ) {
    BatchVector& batches = versioned_batches.back();
    while ( !batches.empty() ) {
      const Batch& batch = batches.back();
      if ( !batch.empty() ) {
        partitioned_hypergraph.uncontract(batch);
      }
      batches.pop_back();
    }
    versioned_batches.pop_back();

    if ( !removed_hyperedges.empty() ) {
      utils::Timer::instance().start_timer(timer_key("restore_parallel_nets"), "Restore Parallel Nets");
      partitioned_hypergraph.restoreSinglePinAndParallelNets(removed_hyperedges.back());
      removed_hyperedges.pop_back();
      utils::Timer::instance().stop_timer(timer_key("restore_parallel_nets"));
    }
  }
  utils::Timer::instance().stop_timer(timer_key("batch_uncontractions"));

  return coarsest_hypergraph;
}

TEST(ANlevel, SimulatesContractionsAndBatchUncontractions) {
  const HypernodeID num_hypernodes = 10000;
  const HypernodeID num_hyperedges = 10000;
  const HypernodeID max_edge_size = 30;
  const HypernodeID num_contractions = 9950;
  const size_t batch_size = 100;
  const bool show_timings = false;
  const bool debug = false;

  if ( debug ) LOG << "Generate Random Hypergraph";
  DynamicHypergraph original_hypergraph = generateRandomHypergraph(num_hypernodes, num_hyperedges, max_edge_size);
  DynamicHypergraph sequential_hg = original_hypergraph.copy(parallel_tag_t());
  DynamicPartitionedHypergraph sequential_phg(4, sequential_hg, parallel_tag_t());
  DynamicHypergraph parallel_hg = original_hypergraph.copy(parallel_tag_t());
  DynamicPartitionedHypergraph parallel_phg(4, parallel_hg, parallel_tag_t());

  if ( debug ) LOG << "Determine random contractions";
  BatchVector contractions = generateRandomContractions(num_hypernodes, num_contractions);

  utils::Timer::instance().clear();

  if ( debug ) LOG << "Simulate n-Level sequentially";
  utils::Timer::instance().start_timer("sequential_n_level", "Sequential n-Level");
  DynamicHypergraph coarsest_sequential_hg = simulateNLevel(sequential_hg, sequential_phg, contractions, 1, false);
  utils::Timer::instance().stop_timer("sequential_n_level");

  if ( debug ) LOG << "Simulate n-Level in parallel";
  utils::Timer::instance().start_timer("parallel_n_level", "Parallel n-Level");
  DynamicHypergraph coarsest_parallel_hg = simulateNLevel(parallel_hg, parallel_phg, contractions, batch_size, true);
  utils::Timer::instance().stop_timer("parallel_n_level");

  if ( debug ) LOG << "Verify equality of hypergraphs";
  verifyEqualityOfHypergraphs(coarsest_sequential_hg, coarsest_parallel_hg);
  verifyEqualityOfHypergraphs(original_hypergraph, sequential_hg);
  verifyEqualityOfHypergraphs(original_hypergraph, parallel_hg);

  if ( debug ) LOG << "Verify gain cache of hypergraphs";
  verifyGainCache(sequential_phg);
  verifyGainCache(parallel_phg);

  if ( debug ) LOG << "Verify number of incident cut hyperedges";
  verifyNumIncidentCutHyperedges(sequential_phg);
  verifyNumIncidentCutHyperedges(parallel_phg);

  if ( show_timings ) {
    LOG << utils::Timer::instance(true);
  }
}

TEST(ANlevel, SimulatesParallelContractionsAndAccessToHypergraph) {
  const HypernodeID num_hypernodes = 10000;
  const HypernodeID num_hyperedges = 10000;
  const HypernodeID max_edge_size = 30;
  const HypernodeID num_contractions = 9950;
  const bool show_timings = false;
  const bool debug = false;

  if ( debug ) LOG << "Generate Random Hypergraph";
  DynamicHypergraph hypergraph = generateRandomHypergraph(num_hypernodes, num_hyperedges, max_edge_size);
  DynamicHypergraph tmp_hypergraph = hypergraph.copy(parallel_tag_t());

  if ( debug ) LOG << "Determine random contractions";
  BatchVector contractions = generateRandomContractions(num_hypernodes, num_contractions, false);

  utils::Timer::instance().clear();

  if ( debug ) LOG << "Perform Parallel Contractions With Parallel Access";
  bool terminate = false;
  utils::Timer::instance().start_timer("contractions_with_access", "Contractions With Access");
  tbb::parallel_invoke([&] {
    while ( !terminate ) {
      // Iterate over all vertices of the hypergraph in parallel
      hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        RatingType rating = 0;
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
          for ( const HypernodeID& pin : hypergraph.pins(he) ) {
            const HyperedgeID node_degree = hypergraph.nodeDegree(pin);
            const HypernodeWeight node_weight = hypergraph.nodeWeight(pin);
            if ( hypergraph.communityID(hn) == hypergraph.communityID(pin) ) {
              rating += static_cast<RatingType>(edge_weight * node_degree) / node_weight;
            }
          }
        }
      });
    }
  }, [&] {
    // Perform contractions in parallel
    tbb::parallel_for(0UL, contractions.back().size(), [&](const size_t i) {
      const Memento& memento = contractions.back()[i];
      hypergraph.registerContraction(memento.u, memento.v);
      hypergraph.contract(memento.v);
    });
    terminate = true;
  });
  utils::Timer::instance().stop_timer("contractions_with_access");

  if ( debug ) LOG << "Perform Parallel Contractions Without Parallel Access";
  utils::Timer::instance().start_timer("contractions_without_access", "Contractions Without Access");
  tbb::parallel_for(0UL, contractions.back().size(), [&](const size_t i) {
    const Memento& memento = contractions.back()[i];
    tmp_hypergraph.registerContraction(memento.u, memento.v);
    tmp_hypergraph.contract(memento.v);
  });
  utils::Timer::instance().stop_timer("contractions_without_access");

  if ( show_timings ) {
    LOG << utils::Timer::instance(true);
  }
}

} // namespace ds
} // namespace mt_kahypar
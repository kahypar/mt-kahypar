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
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
namespace ds {

using Batch = parallel::scalable_vector<Memento>;
using BatchVector = parallel::scalable_vector<Batch>;
using VersionedBatchVector = parallel::scalable_vector<BatchVector>;

void verifyEqualityOfHypergraphs(const DynamicHypergraph& expected_hypergraph,
                                 const DynamicHypergraph& actual_hypergraph) {
  parallel::scalable_vector<HyperedgeID> expected_incident_edges;
  parallel::scalable_vector<HyperedgeID> actual_incident_edges;
  for ( const HypernodeID& hn : expected_hypergraph.nodes() ) {
    ASSERT_TRUE(actual_hypergraph.nodeIsEnabled(hn));
    ASSERT_EQ(expected_hypergraph.nodeWeight(hn), actual_hypergraph.nodeWeight(hn));
    ASSERT_EQ(expected_hypergraph.nodeDegree(hn), actual_hypergraph.nodeDegree(hn));
    for ( const HyperedgeID he : expected_hypergraph.incidentEdges(hn) ) {
      expected_incident_edges.push_back(he);
    }
    for ( const HyperedgeID he : actual_hypergraph.incidentEdges(hn) ) {
      actual_incident_edges.push_back(he);
    }
    std::sort(expected_incident_edges.begin(), expected_incident_edges.end());
    std::sort(actual_incident_edges.begin(), actual_incident_edges.end());
    ASSERT_EQ(expected_incident_edges.size(), actual_incident_edges.size());
    for ( size_t i = 0; i < expected_incident_edges.size(); ++i ) {
      ASSERT_EQ(expected_incident_edges[i], actual_incident_edges[i]);
    }
    expected_incident_edges.clear();
    actual_incident_edges.clear();
  }

  parallel::scalable_vector<HypernodeID> expected_pins;
  parallel::scalable_vector<HypernodeID> actual_pins;
  for ( const HyperedgeID& he : expected_hypergraph.edges() ) {
    for ( const HyperedgeID he : expected_hypergraph.pins(he) ) {
      expected_pins.push_back(he);
    }
    for ( const HyperedgeID he : actual_hypergraph.pins(he) ) {
      actual_pins.push_back(he);
    }
    std::sort(expected_pins.begin(), expected_pins.end());
    std::sort(actual_pins.begin(), actual_pins.end());
    ASSERT_EQ(expected_pins.size(), actual_pins.size());
    for ( size_t i = 0; i < expected_pins.size(); ++i ) {
      ASSERT_EQ(expected_pins[i], actual_pins[i]);
    }
    expected_pins.clear();
    actual_pins.clear();
  }
}

DynamicHypergraph generateRandomHypergraph(const HypernodeID num_hypernodes,
                                           const HyperedgeID num_hyperedges,
                                           const HypernodeID max_edge_size) {
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> hyperedges;
  utils::Randomize& rand = utils::Randomize::instance();
  for ( size_t i = 0; i < num_hyperedges; ++i ) {
    parallel::scalable_vector<HypernodeID> net;
    const size_t edge_size = rand.getRandomInt(2, max_edge_size, sched_getcpu());
    for ( size_t i = 0; i < edge_size; ++i ) {
      const HypernodeID pin = rand.getRandomInt(0, num_hypernodes - 1, sched_getcpu());
      if ( std::find(net.begin(), net.end(), pin) == net.end() ) {
        net.push_back(pin);
      }
    }
    hyperedges.emplace_back(std::move(net));
  }
  return DynamicHypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, num_hypernodes, num_hyperedges, hyperedges);
}

parallel::scalable_vector<Memento> generateRandomContractions(const HypernodeID num_hypernodes,
                                                              const HypernodeID num_contractions) {
 ASSERT(num_contractions < num_hypernodes);
 parallel::scalable_vector<Memento> contractions;
 parallel::scalable_vector<HypernodeID> active_hns(num_hypernodes);
 std::iota(active_hns.begin(), active_hns.end(), 0);
 utils::Randomize& rand = utils::Randomize::instance();
 const int cpu_id = sched_getcpu();
 for ( size_t i = 0; i < num_contractions; ++i ) {
   ASSERT(active_hns.size() >= 2UL);
   int idx_1 = rand.getRandomInt(0, static_cast<int>(active_hns.size() - 1), cpu_id);
   int idx_2 = rand.getRandomInt(0, static_cast<int>(active_hns.size() - 1), cpu_id);
    if ( idx_1 == idx_2 ) {
      idx_2 = (idx_2 + 1) % active_hns.size();
    }
    contractions.push_back(Memento { active_hns[idx_1], active_hns[idx_2] });
    std::swap(active_hns[idx_2], active_hns.back());
    active_hns.pop_back();
 }
 return contractions;
}

DynamicHypergraph simulateNLevel(DynamicHypergraph& hypergraph,
                                 const parallel::scalable_vector<Memento>& contractions,
                                 const size_t batch_size,
                                 const bool parallel) {

  utils::Timer::instance().start_timer("contractions", "Contractions");
  if ( parallel ) {
    tbb::parallel_for(0UL, contractions.size(), [&](const size_t i) {
      const Memento& memento = contractions[i];
      hypergraph.registerContraction(memento.u, memento.v);
      hypergraph.contract(memento.v);
    });
  } else {
    for ( size_t i = 0; i < contractions.size(); ++i ) {
      const Memento& memento = contractions[i];
      hypergraph.registerContraction(memento.u, memento.v);
      hypergraph.contract(memento.v);
    }
  }
  utils::Timer::instance().stop_timer("contractions");

  utils::Timer::instance().start_timer("copy_coarsest_hypergraph", "Copy Coarsest Hypergraph");
  DynamicHypergraph coarsest_hypergraph;
  if ( parallel ) {
    coarsest_hypergraph = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  } else {
    coarsest_hypergraph = hypergraph.copy();
  }
  utils::Timer::instance().stop_timer("copy_coarsest_hypergraph");

  std::string key = "create_batch_uncontraction_hierarchy";
  if ( parallel ) key = "create_parallel_batch_uncontraction_hierarchy";
  utils::Timer::instance().start_timer(key, "Create n-Level Hierarchy");
  const size_t tmp_batch_size = parallel ? batch_size : 1;
  auto versioned_batches = hypergraph.createBatchUncontractionHierarchy(TBBNumaArena::GLOBAL_TASK_GROUP, tmp_batch_size);
  utils::Timer::instance().stop_timer(key);


  utils::Timer::instance().start_timer("batch_uncontractions", "Batch Uncontractions");
  while ( !versioned_batches.empty() ) {
    BatchVector& batches = versioned_batches.back();
    while ( !batches.empty() ) {
      const parallel::scalable_vector<Memento> batch = batches.back();
      hypergraph.uncontract(batch);
      batches.pop_back();
    }
    versioned_batches.pop_back();
  }
  utils::Timer::instance().stop_timer("batch_uncontractions");

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
  DynamicHypergraph sequential_hg = original_hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);
  DynamicHypergraph parallel_hg = original_hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);

  if ( debug ) LOG << "Determine random contractions";
  parallel::scalable_vector<Memento> contractions = generateRandomContractions(num_hypernodes, num_contractions);

  utils::Timer::instance().clear();

  if ( debug ) LOG << "Simulate n-Level sequentially";
  utils::Timer::instance().start_timer("sequential_n_level", "Sequential n-Level");
  DynamicHypergraph coarsest_sequential_hg = simulateNLevel(sequential_hg, contractions, 1, false);
  utils::Timer::instance().stop_timer("sequential_n_level");

  if ( debug ) LOG << "Simulate n-Level in parallel";
  utils::Timer::instance().start_timer("parallel_n_level", "Parallel n-Level");
  DynamicHypergraph coarsest_parallel_hg = simulateNLevel(parallel_hg, contractions, batch_size, true);
  utils::Timer::instance().stop_timer("parallel_n_level");

  if ( debug ) LOG << "Verify equality of hypergraphs";
  verifyEqualityOfHypergraphs(coarsest_sequential_hg, coarsest_parallel_hg);
  verifyEqualityOfHypergraphs(original_hypergraph, sequential_hg);
  verifyEqualityOfHypergraphs(original_hypergraph, parallel_hg);

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
  DynamicHypergraph tmp_hypergraph = hypergraph.copy(TBBNumaArena::GLOBAL_TASK_GROUP);

  if ( debug ) LOG << "Determine random contractions";
  parallel::scalable_vector<Memento> contractions = generateRandomContractions(num_hypernodes, num_contractions);

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
    tbb::parallel_for(0UL, contractions.size(), [&](const size_t i) {
      const Memento& memento = contractions[i];
      hypergraph.registerContraction(memento.u, memento.v);
      hypergraph.contract(memento.v);
    });
    terminate = true;
  });
  utils::Timer::instance().stop_timer("contractions_with_access");

  if ( debug ) LOG << "Perform Parallel Contractions Without Parallel Access";
  utils::Timer::instance().start_timer("contractions_without_access", "Contractions Without Access");
  tbb::parallel_for(0UL, contractions.size(), [&](const size_t i) {
    const Memento& memento = contractions[i];
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
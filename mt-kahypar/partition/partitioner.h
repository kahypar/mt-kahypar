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

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include "kahypar/meta/policy_registry.h"

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/single_node_hyperedge_remover.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {
namespace partition {

class Partitioner {
 private:
  static constexpr bool debug = true;

 public:
  Partitioner() :
    _single_node_he_remover() { }

  Partitioner(const Partitioner&) = delete;
  Partitioner& operator= (const Partitioner&) = delete;

  Partitioner(Partitioner&&) = delete;
  Partitioner& operator= (Partitioner&&) = delete;

  ~Partitioner() = default;

  inline void partition(Hypergraph& hypergraph, Context& context);

 private:

  static inline void setupContext(const Hypergraph& hypergraph, Context& context);

  static inline void configurePreprocessing(const Hypergraph& hypergraph, Context& context);

  inline void sanitize(Hypergraph& hypergraph, const Context& context);

  inline void preprocess(Hypergraph& hypergraph, const Context& context);

  inline void redistribution(Hypergraph& hypergraph, const Context& context);

  inline void postprocess(Hypergraph& hypergraph);

  SingleNodeHyperedgeRemover _single_node_he_remover;
};

inline void Partitioner::setupContext(const Hypergraph& hypergraph, Context& context) {
  unused(hypergraph);
  unused(context);
}

inline void Partitioner::configurePreprocessing(const Hypergraph& hypergraph, Context& context) {
  unused(hypergraph);
  unused(context);
}

inline void Partitioner::sanitize(Hypergraph& hypergraph, const Context& context) {
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  const auto result = _single_node_he_remover.removeSingleNodeHyperedges(hypergraph);
  if (context.partition.verbose_output && result.num_removed_single_node_hes > 0) {
    LOG << "Performing single-node HE removal:";
    LOG << "\033[1m\033[31m" << " # removed hyperedges with |e|=1 = "
        << result.num_removed_single_node_hes
        << "\033[0m";
    LOG << "\033[1m\033[31m" << " ===>" << result.num_unconnected_hns
        << "unconnected HNs could have been removed" << "\033[0m";
    io::printStripe();
  }
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("single_node_hyperedge_removal", "Single Node Hyperedge Removal",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 0, std::chrono::duration<double>(end - start).count());
}

inline void Partitioner::preprocess(Hypergraph& hypergraph, const Context& context) {
  HighResClockTimepoint global_start = std::chrono::high_resolution_clock::now();
  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  std::vector<PartitionID> communities;
  io::readPartitionFile(context.partition.graph_community_filename, communities);
  ASSERT(communities.size() == hypergraph.initialNumNodes());
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("read_community_file", "Read Community File",
    "community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING, 0, std::chrono::duration<double>(end - start).count());

  // Stream community ids into hypergraph
  start = std::chrono::high_resolution_clock::now();
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
    [&](const tbb::blocked_range<HypernodeID>& range) {
    for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
      hypergraph.streamCommunityID(hypergraph.globalNodeID(hn), communities[hn]);
    }
  });  
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("stream_community_ids", "Stream Community IDs",
    "community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

  // Initialize Communities
  start = std::chrono::high_resolution_clock::now();
  hypergraph.initializeCommunities();
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_communities", "Initialize Communities",
    "community_detection", mt_kahypar::utils::Timer::Type::PREPROCESSING, 2, std::chrono::duration<double>(end - start).count());
  HighResClockTimepoint global_end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("community_detection", "Community Detection",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(global_end - global_start).count());

  if ( context.shared_memory.use_community_redistribution ) {
    // Redistribute Hypergraph based on communities
    start = std::chrono::high_resolution_clock::now();
    redistribution(hypergraph, context);
    end = std::chrono::high_resolution_clock::now();
    mt_kahypar::utils::Timer::instance().add_timing("redistribution", "Redistribution",
      "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 2, std::chrono::duration<double>(end - start).count());
  }
}

inline void Partitioner::redistribution(Hypergraph& hypergraph, const Context& context) {
  std::unique_ptr<preprocessing::IRedistribution> redistributor =
    RedistributionFactory::getInstance().createObject(
      context.shared_memory.assignment_strategy, hypergraph, context);

  DBG << "Remote Pin Count Before Redistribution" << metrics::remotePinCount(hypergraph);
  hypergraph = redistributor->redistribute();
  DBG << "Remote Pin Count After Redistribution" << metrics::remotePinCount(hypergraph);

}

inline void Partitioner::postprocess(Hypergraph& hypergraph) {
  _single_node_he_remover.restoreSingleNodeHyperedges(hypergraph);
}

inline void Partitioner::partition(Hypergraph& hypergraph, Context& context) {
  configurePreprocessing(hypergraph, context);

  setupContext(hypergraph, context);
  io::printInputInformation(context, hypergraph);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  preprocess(hypergraph, context);
  sanitize(hypergraph, context);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("preprocessing", "Preprocessing",
    "", mt_kahypar::utils::Timer::Type::PREPROCESSING, 1, std::chrono::duration<double>(end - start).count());

  // partition
  start = std::chrono::high_resolution_clock::now();
  hypergraph.initializeCommunityHyperedges();  
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("initialize_community_hyperedges", "Initialize Community Hyperedges",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 3, std::chrono::duration<double>(end - start).count());

  start = std::chrono::high_resolution_clock::now();
  hypergraph.resetCommunityHyperedges( {} );
  end = std::chrono::high_resolution_clock::now();
  mt_kahypar::utils::Timer::instance().add_timing("reset_community_hyperedges", "Reset Community Hyperedges",
    "preprocessing", mt_kahypar::utils::Timer::Type::PREPROCESSING, 4, std::chrono::duration<double>(end - start).count());


  postprocess(hypergraph);
}


} // namespace partition
} // namespace mt_kahypar
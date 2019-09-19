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

#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/preprocessing/single_node_hyperedge_remover.h"

namespace mt_kahypar {
namespace partition {

class Partitioner {
 private:
  static constexpr bool debug = false;

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
}

inline void Partitioner::preprocess(Hypergraph& hypergraph, const Context& context) {
  unused(hypergraph);
  std::vector<PartitionID> communities;
  io::readPartitionFile(context.partition.graph_community_filename, communities);
  ASSERT(communities.size() == hypergraph.initialNumNodes());
  
  // Stream community ids into hypergraph
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
    [&](const tbb::blocked_range<HypernodeID>& range) {
    for ( HypernodeID hn = range.begin(); hn < range.end(); ++hn ) {
      hypergraph.streamCommunityID(hypergraph.globalNodeID(hn), communities[hn]);
    }
  });
}

inline void Partitioner::postprocess(Hypergraph& hypergraph) {
  _single_node_he_remover.restoreSingleNodeHyperedges(hypergraph);
}

inline void Partitioner::partition(Hypergraph& hypergraph, Context& context) {
  configurePreprocessing(hypergraph, context);

  setupContext(hypergraph, context);
  io::printInputInformation(context, hypergraph);

  sanitize(hypergraph, context);

  preprocess(hypergraph, context);
  // partition
  postprocess(hypergraph);
}


} // namespace partition
} // namespace mt_kahypar
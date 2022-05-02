/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "tbb/parallel_sort.h"

#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/registries/register_memory_pool.h"
#include "mt-kahypar/io/hypergraph_io.h"

using namespace mt_kahypar;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  Context context;
  mt_kahypar::processCommandLineInput(context, argc, argv);
  context.partition.num_vcycles = 1;

  mt_kahypar::TBBInitializer::instance(context.shared_memory.num_threads);

  // Read Hypergraph
  Hypergraph hg = mt_kahypar::io::readInputFile(
    context.partition.graph_filename, context.partition.file_format, true, false);
  ALWAYS_ASSERT(hg.totalWeight() == static_cast<HypernodeWeight>(hg.initialNumNodes()));

  mt_kahypar::register_memory_pool(hg, context);

  context.partition.large_hyperedge_size_threshold = std::max(hg.initialNumNodes() *
                                                              context.partition.large_hyperedge_size_threshold_factor, 100.0);
  context.sanityCheck();
  context.setupPartWeights(hg.totalWeight());
  context.setupContractionLimit(hg.totalWeight());
  context.setupThreadsPerFlowSearch();

  parallel::scalable_vector<std::pair<HypernodeID, HyperedgeID>> hns_with_degree;
  hns_with_degree.resize(hg.initialNumNodes());
  tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID hn) {
    hns_with_degree[hn] = {hn, hg.nodeDegree(hn)};
  });
  tbb::parallel_sort(hns_with_degree.begin(), hns_with_degree.end(),
    [&](const auto& left, const auto& right) {
      return left.second < right.second;
    }
  );

  PartitionedHypergraph phg(2, hg);
  tbb::parallel_for(ID(0), hg.initialNumNodes(), [&](HypernodeID hn) {
    PartitionID part = (hn <= hg.initialNumNodes() / 2) ? 0 : 1;
    phg.setOnlyNodePart(hns_with_degree[hn].first, part);
  });
  phg.initializePartition();

  HyperedgeWeight cut = metrics::hyperedgeCut(phg);

  LOG << "";
  LOG << "initial_cut=" << cut;

  multilevel::partitionVCycle(hg, phg, context);

  cut = metrics::hyperedgeCut(phg);

  std::string graph_name = context.partition.graph_filename.substr(
    context.partition.graph_filename.find_last_of("/") + 1);
  std::cout  << "RESULT graph=" << graph_name
             << " HNs=" << hg.initialNumNodes()
             << " HEs=" << hg.initialNumEdges()
             << " pins=" << hg.initialNumPins()
             << " density=" << static_cast<double>(hg.initialNumEdges()) / hg.initialNumNodes()
             << " medianDegree=" << hns_with_degree[hg.initialNumNodes() / 2].second
             << " cut=" << cut
             << std::endl;

  return 0;
}

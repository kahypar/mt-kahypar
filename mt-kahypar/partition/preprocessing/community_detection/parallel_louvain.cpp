/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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


#include "parallel_louvain.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar::community_detection {

  template<typename Hypergraph>
  std::vector<std::pair<ds::Clustering, double>> local_moving_contract_recurse(Graph<Hypergraph>& fine_graph,
                                                                               ParallelLocalMovingModularity<Hypergraph>& mlv,
                                                                               const Context& context) {
    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
    timer.start_timer("local_moving", "Local Moving");
    ds::Clustering own_communities(fine_graph.numNodes());
    bool communities_changed = mlv.localMoving(fine_graph, own_communities);
    timer.stop_timer("local_moving");
    std::vector<std::pair<ds::Clustering, double>> result;

    if (communities_changed) {
      timer.start_timer("contraction_cd", "Contraction");
      // Contract Communities
      Graph<Hypergraph> coarse_graph = fine_graph.contract(own_communities, context.preprocessing.community_detection.low_memory_contraction);
      ASSERT(coarse_graph.totalVolume() == fine_graph.totalVolume());
      timer.stop_timer("contraction_cd");

      double new_modularity = 0;
      double factor = 1 / coarse_graph.totalVolume();
      for (NodeID node: coarse_graph.nodes()) {
        double contribution = coarse_graph.nodeVolume(node);
        for (const Arc& arc : coarse_graph.arcsOf(node)) {
          contribution -= arc.weight;  // only count internal edges
        }
        contribution -= factor * coarse_graph.nodeVolume(node) * coarse_graph.nodeVolume(node);
        new_modularity += factor * contribution;
      }
      result.emplace_back(own_communities, new_modularity);

      // Recurse on contracted graph
      auto coarse_communities = local_moving_contract_recurse(coarse_graph, mlv, context);

      timer.start_timer("project", "Project");
      // Prolong Clustering
      for (const auto& [comm, modularity]: coarse_communities) {
        ds::Clustering communities(own_communities);  // yes, this is an intentional copy
        tbb::parallel_for(UL(0), fine_graph.numNodes(), [&](const NodeID u) {
          ASSERT(communities[u] < static_cast<PartitionID>(comm.size()));
          communities[u] = comm[communities[u]];
        });
        result.emplace_back(std::move(communities), modularity);
      }
      timer.stop_timer("project");
    }

    return result;
  }

  template<typename Hypergraph>
  std::vector<std::pair<ds::Clustering, double>> run_parallel_louvain(Graph<Hypergraph>& graph,
                                      const Context& context,
                                      bool disable_randomization) {
    ParallelLocalMovingModularity<Hypergraph> mlv(context, graph.numNodes(), disable_randomization);
    auto result = local_moving_contract_recurse(graph, mlv, context);
    return result;
  }

  namespace {
  #define LOCAL_MOVING(X) std::vector<std::pair<ds::Clustering, double>> local_moving_contract_recurse(Graph<X>&, ParallelLocalMovingModularity<X>&, const Context&)
  #define PARALLEL_LOUVAIN(X) std::vector<std::pair<ds::Clustering, double>> run_parallel_louvain(Graph<X>&, const Context&, bool)
  }

  INSTANTIATE_FUNC_WITH_HYPERGRAPHS(LOCAL_MOVING)
  INSTANTIATE_FUNC_WITH_HYPERGRAPHS(PARALLEL_LOUVAIN)
}
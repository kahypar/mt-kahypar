/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <algorithm>

#include "mt-kahypar/partition/refinement/rebalancing/kaminpar_rebalancer.h"

#include <common/datastructures/static_array.h>
#include <common/parallel/algorithm.h>
#include <common/timer.h>
#include <kaminpar/context.h>
#include <kaminpar/datastructures/graph.h>
#include <kaminpar/datastructures/partitioned_graph.h>
#include <kaminpar/kaminpar.h>
#include <kaminpar/metrics.h>

#include "tbb/parallel_for.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/partition/refinement/rebalancing/call_kaminpar.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/randomize.h"


namespace mt_kahypar {
  template <typename TypeTraits, typename GainTypes>
  bool KaminparRebalancer<TypeTraits, GainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                            const vec<HypernodeID>&,
                                                            Metrics& best_metrics,
                                                            double) {
    using namespace kaminpar;
    using namespace kaminpar::shm;

    PartitionedHypergraph& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    if (_max_weight == 0) {
      _max_weight = _context.partition.max_part_weights[0] * _context.partition.k;
    }
    const double epsilon = 1.0 * _max_weight / phg.totalWeight() - 1.0;

    // If partition is imbalanced, rebalancer is activated
    bool improvement = false;
    if ( !isBalanced(phg) ) {
      DISABLE_TIMERS();
      _gain_cache.reset();

      DBG << "[MtKaHyPar] Metrics *before* calling KaminparRebalancer: cut="
          << metrics::quality(phg, _context)
          << " imbalance=" << metrics::imbalance(phg, _context);
      DBG << "[MtKaHyPar] Number of removed nodes: " << phg.numRemovedHypernodes();

      const NodeID n =
          phg.initialNumNodes() - phg.numRemovedHypernodes();
      const EdgeID m = phg.initialNumEdges();

      StaticArray<NodeID> dense(phg.initialNumNodes() + 1);
      phg.doParallelForAllNodes([&](const HypernodeID u) {
          KASSERT(u + 1 < dense.size());
          dense[u + 1] = 1;
      });
      kaminpar::parallel::prefix_sum(dense.begin(), dense.end(), dense.begin());
      KASSERT(dense.front() == 0u);
      KASSERT(dense.back() == n);

      StaticArray<EdgeID> xadj(n + 1);
      StaticArray<NodeID> adjncy(m);
      StaticArray<NodeWeight> vwgt(n);
      StaticArray<EdgeWeight> adjwgt(m);
      StaticArray<BlockID> part(n);

      phg.doParallelForAllNodes([&](const HypernodeID u) {
        const NodeID du = dense[u];
        xadj[du + 1] = phg.nodeDegree(u);

        const HypernodeWeight wgt = phg.nodeWeight(u);
        vwgt[du] = wgt;

        const PartitionID block = phg.partID(u);
        part[du] = block;
      });

      kaminpar::parallel::prefix_sum(xadj.begin(), xadj.end(), xadj.begin());

      phg.doParallelForAllNodes([&](const HypernodeID u) {
        const NodeID du = dense[u];
        HypernodeID offset = 0;

        for (const HyperedgeID e : phg.incidentEdges(u)) {
          const HypernodeID v = phg.edgeTarget(e);
          const NodeID dv = dense[v];

          const std::size_t index = xadj[du] + offset;
          adjncy[index] = dv;
          adjwgt[index] = phg.edgeWeight(e);

          ++offset;
        }
      });

      kaminpar::shm::Graph graph(std::move(xadj), std::move(adjncy),
                                std::move(vwgt), std::move(adjwgt), false);
      kaminpar::shm::PartitionedGraph p_graph(graph, phg.k(),
                                              std::move(part));

      kaminpar::shm::Context ctx = kaminpar::shm::create_default_context();
      ctx.partition.k = phg.k();
      ctx.partition.epsilon = epsilon;
      ctx.setup(graph);

      DBG << "[KaMinPar] Metrics *before* calling KaminparRebalancer: cut="
          << kaminpar::shm::metrics::edge_cut(p_graph)
          << " imbalance=" << kaminpar::shm::metrics::imbalance(p_graph)
          << " epsilon=" << ctx.partition.epsilon
          << " max_block_weight=" << ctx.partition.block_weights.max(0)
          << " max_block_weight'=" << _context.partition.max_part_weights[0];

      kaminpar::call_balancer(p_graph, ctx);

      DBG << "[KaMinPar] Metrics *after* calling KaminparRebalancer: cut="
          << kaminpar::shm::metrics::edge_cut(p_graph)
          << " imbalance=" << kaminpar::shm::metrics::imbalance(p_graph);

      phg.doParallelForAllNodes([&](const HypernodeID u) {
        const NodeID du = dense[u];

        const PartitionID move_from = phg.partID(u);
        const PartitionID move_to = p_graph.block(du);

        if (move_from != move_to) {
          phg.changeNodePart(u, move_from, move_to);
        }
      });

      DBG << "[MtKaHyPar] Metrics *after* calling KaminparRebalancer: cut="
          << metrics::quality(phg, _context)
          << " imbalance=" << metrics::imbalance(phg, _context);

      HyperedgeWeight new_quality = metrics::quality(phg, _context);
      improvement = new_quality < best_metrics.quality;
      best_metrics.quality = new_quality;
      best_metrics.imbalance = metrics::imbalance(phg, _context);
    }
    _max_weight = 0;
    return improvement;;
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  namespace {
  #define KAMINPAR_REBALANCER(X, Y) KaminparRebalancer<X, Y>
  }

  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(KAMINPAR_REBALANCER)
}
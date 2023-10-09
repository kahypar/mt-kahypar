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


#pragma once

#include <algorithm>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar-resources/meta/policy_registry.h"
#include "kahypar-resources/meta/typelist.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

class AlwaysAcceptPolicy final : public kahypar::meta::PolicyBase {
 public:
  AlwaysAcceptPolicy() { }

  explicit AlwaysAcceptPolicy(const HypernodeID) { }

  template<typename Hypergraph>
  void initialize(const Hypergraph&, const Context&) { }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool acceptContraction(const Hypergraph&, const Context&, HypernodeID, HypernodeID) const {
    return true;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  double similarityPenalty(const Hypergraph&, const Context&, HypernodeID, HypernodeID) const {
    return 1.0;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  double weightRatioForNode(const Hypergraph& hypergraph, const HypernodeID u) const {
    return 1.0;
  }
};


class PreserveRebalancingNodesPolicy final : public kahypar::meta::PolicyBase {
  static constexpr bool debug = false;

 public:
  PreserveRebalancingNodesPolicy(): PreserveRebalancingNodesPolicy(0) {}

  explicit PreserveRebalancingNodesPolicy(const HypernodeID num_nodes):
    _incident_weight(num_nodes, 0), _acceptance_limit(num_nodes, 0) {}

  PreserveRebalancingNodesPolicy(const PreserveRebalancingNodesPolicy&) = delete;
  PreserveRebalancingNodesPolicy(PreserveRebalancingNodesPolicy&&) = delete;
  PreserveRebalancingNodesPolicy & operator= (const PreserveRebalancingNodesPolicy &) = delete;
  PreserveRebalancingNodesPolicy & operator= (PreserveRebalancingNodesPolicy &&) = delete;

  template<typename Hypergraph>
  void initialize(const Hypergraph& hypergraph, const Context& context) {
    ASSERT(_incident_weight.size() >= hypergraph.initialNumNodes()
           && _acceptance_limit.size() >= hypergraph.initialNumNodes());

    struct NeighborData {
      HypernodeID hn;
      float edge_weight_contribution;
      HypernodeWeight node_weight;
    };

    auto scaled_edge_weight = [&](const HyperedgeID he) {
      if constexpr (Hypergraph::is_graph) {
        return hypergraph.edgeWeight(he);
      } else {
        return static_cast<double>(hypergraph.edgeWeight(he)) /
          (hypergraph.edgeSize(he) + context.coarsening.rating.incident_weight_scaling_constant);
      }
    };

    // compute incident weights
    hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      // TODO(maas): save the total incident weight in the hypergraph data structure?
      double incident_weight_sum = 0;
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        incident_weight_sum += scaled_edge_weight(he);
      }
      _incident_weight[hn] = incident_weight_sum;
    });
    if constexpr (Hypergraph::is_graph) {
      // TODO: We are ignoring edges between neighbors here - the result is thus only approximate.
      // This could be acceptable, though
      const HypernodeWeight max_summed_weight = std::ceil(context.coarsening.rating.preserve_nodes_relative_weight_limit
                                                          * hypergraph.totalWeight());
      tbb::enumerable_thread_specific<parallel::scalable_vector<NeighborData>> local_neighbor_list;
      hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        auto& neighbor_list = local_neighbor_list.local();
        neighbor_list.clear();

        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          HypernodeID v = hypergraph.edgeTarget(he);
          // TODO: we could filter the nodes that can not decrease the value
          neighbor_list.push_back(NeighborData{v, _incident_weight[v] - 2 * scaled_edge_weight(he), hypergraph.nodeWeight(v)});
        }
        std::sort(neighbor_list.begin(), neighbor_list.end(), [](const NeighborData& lhs, const NeighborData& rhs) {
          if (lhs.node_weight == 0) {
            return false;
          } else if (rhs.node_weight == 0) {
            return true;
          }
          return static_cast<double>(lhs.edge_weight_contribution) / lhs.node_weight
            < static_cast<double>(rhs.edge_weight_contribution) / rhs.node_weight;
        });

        const HypernodeID max_iterations = context.coarsening.rating.max_considered_neighbors;
        double summed_contribution = _incident_weight[hn];
        HypernodeWeight summed_weight = std::max(hypergraph.nodeWeight(hn), 1);
        double current_min = summed_contribution / summed_weight;
        for (size_t i = 0; i < neighbor_list.size() && (max_iterations == 0 || i < max_iterations); ++i) {
          const NeighborData& neighbor = neighbor_list[i];
          summed_contribution += neighbor.edge_weight_contribution;
          summed_weight += neighbor.node_weight;
          if (summed_contribution / summed_weight <= current_min && summed_weight <= max_summed_weight) {
            current_min = summed_contribution / summed_weight;
          } else {
            break;
          }
        }
        _acceptance_limit[hn] = std::min(
          context.coarsening.rating.preserve_nodes_scaling_factor * current_min,
          context.coarsening.rating.acceptance_limit_bound * _incident_weight[hn] / std::max(hypergraph.nodeWeight(hn), 1));
        DBG << V(hn) << V(_acceptance_limit[hn]) << V(_incident_weight[hn])
            << V(hypergraph.nodeWeight(hn)) << V(hypergraph.nodeDegree(hn));
      });
    } else {
      ERR("not supported");
    }
  }

  // this function decides if contracting v onto u is allowed
  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool acceptContraction(const Hypergraph& hypergraph,
                                                            const Context&,
                                                            const HypernodeID u,
                                                            const HypernodeID v) const {
    double ratio_u = weightRatioForNode(hypergraph, u);
    double ratio_v = weightRatioForNode(hypergraph, v);
    if (ratio_u >= ratio_v) {
      return ratio_v >= _acceptance_limit[u] || hypergraph.nodeWeight(v) == 0;
    } else {
      return ratio_u >= _acceptance_limit[v] || hypergraph.nodeWeight(u) == 0;
    }
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE double similarityPenalty(const Hypergraph& hypergraph,
                                                              const Context& context,
                                                              const HypernodeID u,
                                                              const HypernodeID v) const {
    if (context.coarsening.rating.use_similarity_penalty) {
      double ratio_u = weightRatioForNode(hypergraph, u);
      double ratio_v = weightRatioForNode(hypergraph, v);
      return std::max(ratio_u / ratio_v, ratio_v / ratio_u);
    }
    return 1.0;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE double weightRatioForNode(const Hypergraph& hypergraph,
                                                               const HypernodeID u) const {
    return _incident_weight[u] / std::max(hypergraph.nodeWeight(u), 1);
  }

 private:
  parallel::scalable_vector<float> _incident_weight;   // TODO: use ints for graphs??
  parallel::scalable_vector<float> _acceptance_limit;
};


using DegreeSimilarityPolicies = kahypar::meta::Typelist<PreserveRebalancingNodesPolicy>;

}  // namespace mt_kahypar

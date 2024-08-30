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
  double weightRatioForNode(const Hypergraph& hypergraph, const HypernodeID u) const {
    return 1.0;
  }
};

namespace {
  struct NeighborData {
    float edge_weight_contribution;
    HypernodeWeight node_weight;
  };

  class GroupedIncidenceData {
   public:
    static constexpr size_t MAX_NUM_GROUPS = 16;
    static constexpr double GROUP_FACTOR = 1.5;

    using GroupedData = std::array<NeighborData, MAX_NUM_GROUPS>;

    GroupedIncidenceData() {
      reset();
    }

    void insert(float edge_contribution, HypernodeWeight weight) {
      ASSERT(weight > 0);
      if (edge_contribution <= 0) {
        _data[0].edge_weight_contribution += edge_contribution;
        _data[0].node_weight += weight;
        return;
      } else if (weight == 0) {
        return;
      }
      const double ratio = static_cast<double>(edge_contribution) / static_cast<double>(weight);

      size_t i = 1;
      for (; i < _data.size(); ++i) {
        auto [current_contribution, current_weight] = _data[i];
        if (current_weight == 0) {
          _data[i] = NeighborData{edge_contribution, weight};
          return;
        }

        double current_ratio = static_cast<double>(current_contribution) / static_cast<double>(current_weight);
        if (ratio <= GROUP_FACTOR * current_ratio) {
          if (current_ratio <= GROUP_FACTOR * ratio) {
            // add the values to the fitting group
            _data[i].edge_weight_contribution += edge_contribution;
            _data[i].node_weight += weight;
          } else {
            // insert a new group and shift all following groups
            NeighborData tmp = NeighborData{edge_contribution, weight};
            for (size_t j = i; tmp.node_weight > 0 && j < _data.size(); ++j) {
              std::swap(tmp, _data[j]);
            }
          }
          return;
        }
      }
    }

    const GroupedData& inner() const {
      return _data;
    }

    void reset() {
      _data.fill(NeighborData{0, 0});
    }

   private:
    GroupedData _data;
  };
} // namespace


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

      // Step 1: Collect contributed edge weights and node weights of neighbors in into sorted aggregates
      // (effectively a semi-sorting)
      const HypernodeWeight max_summed_weight = std::ceil(context.coarsening.rating.preserve_nodes_relative_weight_limit
                                                          * hypergraph.totalWeight());
      hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
        GroupedIncidenceData incidence_data;
        const double ratio_of_u = _incident_weight[hn] / std::max(hypergraph.nodeWeight(hn), 1);
        // TODO: this needs to be implemented differently for hypergraphs
        // TODO: don't consider all neighbors for nodes with very high degree?
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          HypernodeID v = hypergraph.edgeTarget(he);
          float edge_contribution = _incident_weight[v] - 2 * scaled_edge_weight(he);
          HypernodeWeight weight = hypergraph.nodeWeight(v);
          if (weight == 0 || edge_contribution / weight < ratio_of_u) {
            incidence_data.insert(edge_contribution, weight);
          }
        }

        // Step 2: Iterate through aggregated neighbor values in sorted order and determine minimum
        const auto& list = incidence_data.inner();
        double summed_contribution = _incident_weight[hn];
        HypernodeWeight summed_weight = std::max(hypergraph.nodeWeight(hn), 1);
        double min_value = summed_contribution / summed_weight;
        for (size_t i = 0; i < list.size() && summed_weight <= max_summed_weight; ++i) {
          const NeighborData& neighbor = list[i];
          if (summed_weight + neighbor.node_weight > max_summed_weight) {
            double fraction_of_last = static_cast<double>(max_summed_weight - summed_weight) / neighbor.node_weight;
            summed_contribution += fraction_of_last * neighbor.edge_weight_contribution;
            summed_weight = max_summed_weight;
          } else {
            summed_contribution += neighbor.edge_weight_contribution;
            summed_weight += neighbor.node_weight;
          }
          if (summed_contribution / summed_weight <= min_value) {
            min_value = summed_contribution / summed_weight;
          } else {
            break;
          }
        }

        // Step 3: Compute acceptance limit of v from minimum
        _acceptance_limit[hn] = std::min(
          context.coarsening.rating.preserve_nodes_scaling_factor * min_value,
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

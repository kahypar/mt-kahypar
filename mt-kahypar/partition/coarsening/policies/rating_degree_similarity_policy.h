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
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

class AlwaysAcceptPolicy final : public kahypar::meta::PolicyBase {
 public:
  AlwaysAcceptPolicy() { }

  explicit AlwaysAcceptPolicy(const HypernodeID) { }

  template<typename Hypergraph>
  void initialize(const Hypergraph&, const Context&, utils::Timer&) { }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool acceptContraction(const Hypergraph&, const Context&, HypernodeID, HypernodeID) const {
    return true;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void accumulate(const Context&, EdgeMetadata&, const EdgeMetadata, const HyperedgeWeight) const { }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool acceptEdgeContraction(const Hypergraph&,
                                                                const Context&,
                                                                const double,
                                                                const HyperedgeWeight,
                                                                const EdgeMetadata) const {
    return true;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeWeight scaledRating(const Hypergraph&,
                                                                  const Context&,
                                                                  const double,
                                                                  const HyperedgeWeight rating,
                                                                  const EdgeMetadata) const {
    return rating;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  double weightRatioForNode(const Hypergraph&, HypernodeID) const {
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
      ASSERT(weight >= 0);
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
  using IncidenceMap = ds::SparseMap<HypernodeID, float>;  // this is prototypical and will almost certainly be removed
  static constexpr bool debug = false;

 public:
  explicit PreserveRebalancingNodesPolicy():
   _incident_weight(), _acceptance_limit(), _local_incidence_map(0) {}

  explicit PreserveRebalancingNodesPolicy(const HypernodeID num_nodes):
    _incident_weight(num_nodes, 0), _acceptance_limit(num_nodes, 0), _local_incidence_map(num_nodes) {}

  PreserveRebalancingNodesPolicy(const PreserveRebalancingNodesPolicy&) = delete;
  PreserveRebalancingNodesPolicy(PreserveRebalancingNodesPolicy&&) = delete;
  PreserveRebalancingNodesPolicy & operator= (const PreserveRebalancingNodesPolicy &) = delete;
  PreserveRebalancingNodesPolicy & operator= (PreserveRebalancingNodesPolicy &&) = delete;

  template<typename Hypergraph>
  void initialize(const Hypergraph& hypergraph, const Context& context, utils::Timer& timer) {
    ASSERT(_incident_weight.size() >= hypergraph.initialNumNodes()
           && _acceptance_limit.size() >= hypergraph.initialNumNodes());

    auto edge_weight_scaling = [&](const HyperedgeID he) {
      if constexpr (Hypergraph::is_graph) {
        return 1.0;
      } else if (hypergraph.edgeSize(he) <= context.coarsening.rating.incident_weight_scaling_constant) {
        return 1.0;
      } else {
        return context.coarsening.rating.incident_weight_scaling_constant / static_cast<double>(hypergraph.edgeSize(he));
      }
    };

    timer.start_timer("compute_incident_weight", "Compute Incident Weight");
    // compute incident weights
    hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      // TODO(maas): save the total incident weight in the hypergraph data structure?
      double incident_weight_sum = 0;
      for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
        incident_weight_sum += edge_weight_scaling(he) * hypergraph.edgeWeight(he);
      }
      _incident_weight[hn] = incident_weight_sum;
    });
    timer.stop_timer("compute_incident_weight");

    timer.start_timer("compute_similarity_metric", "Compute Similarity Metric");
    // TODO: We are ignoring edges between neighbors here - the result is thus only approximate.
    // This could be acceptable, though
    const HypernodeWeight max_summed_weight = std::ceil(context.coarsening.rating.preserve_nodes_relative_weight_limit
                                                        * hypergraph.totalWeight());
    hypergraph.doParallelForAllNodes([&](const HypernodeID hn) {
      GroupedIncidenceData incidence_data;
      const double ratio_of_u = _incident_weight[hn] / std::max(hypergraph.nodeWeight(hn), 1);

      // Step 1: Collect contributed edge weights and node weights of neighbors in into sorted aggregates
      // (effectively a semi-sorting)
      // TODO: should this rather be relative to the maximum cluster weight?
      if constexpr (Hypergraph::is_graph) {
        size_t num_accesses = 0;
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          HypernodeID v = hypergraph.edgeTarget(he);
          float edge_contribution = _incident_weight[v] - 2 * hypergraph.edgeWeight(he);
          HypernodeWeight weight = hypergraph.nodeWeight(v);
          if (weight == 0 || edge_contribution / weight < ratio_of_u) {
            incidence_data.insert(edge_contribution, weight);
          }
          // TODO: should we scale the threshold down in comparison to the rater?
          ++num_accesses;
          if (num_accesses + 1 > context.coarsening.vertex_degree_sampling_threshold) {
            break;
          }
        }
      } else {
        // this is probably quite slow and will be replaced with a bloom-filter based approach
        size_t num_accesses = 0;
        IncidenceMap& incidence_map = _local_incidence_map.local();
        incidence_map.clear();
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          HypernodeID edge_size = hypergraph.edgeSize(he);
          if (edge_size < context.partition.ignore_hyperedge_size_threshold) {
            if (num_accesses + edge_size > context.coarsening.vertex_degree_sampling_threshold) {
              break;
            }
            for (const HypernodeID& pin: hypergraph.pins(he)) {
              if (pin != hn) {
                incidence_map[pin] += edge_weight_scaling(he) * static_cast<double>(hypergraph.edgeWeight(he)) / (edge_size - 1);
                ++num_accesses;
              }
            }
          }
        }

        for (const auto& [neighbor, connectivity]: incidence_map) {
          float edge_contribution = _incident_weight[neighbor] - 2 * connectivity;
          HypernodeWeight weight = hypergraph.nodeWeight(neighbor);
          if (weight == 0 || edge_contribution / weight < ratio_of_u) {
            incidence_data.insert(edge_contribution, weight);
          }
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
    timer.stop_timer("compute_similarity_metric");
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

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void accumulate(const Context&, EdgeMetadata&, const EdgeMetadata, const HyperedgeWeight) const { }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool acceptEdgeContraction(const Hypergraph&,
                                                                const Context&,
                                                                const double,
                                                                const HyperedgeWeight,
                                                                const EdgeMetadata) const {
    return true;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeWeight scaledRating(const Hypergraph&,
                                                                  const Context&,
                                                                  const double,
                                                                  const HyperedgeWeight rating,
                                                                  const EdgeMetadata) const {
    return rating;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE double weightRatioForNode(const Hypergraph& hypergraph,
                                                               const HypernodeID u) const {
    return _incident_weight[u] / std::max(hypergraph.nodeWeight(u), 1);
  }

 private:
  const Context* _context;  // TODO: currently must be a pointer so we can default-initialize..
  // ! incident weight (scaled with hyperedge size) for all nodes
  parallel::scalable_vector<float> _incident_weight;
  // ! pre-computed metric which is used to determine whether a contraction is accepted
  parallel::scalable_vector<float> _acceptance_limit;
  // ! Tracks connectivity to all neighbors in case of hypergraphs
  tbb::enumerable_thread_specific<IncidenceMap> _local_incidence_map;
};



class GuidedCoarseningPolicy final : public kahypar::meta::PolicyBase {
  static constexpr bool debug = false;

 public:
  explicit GuidedCoarseningPolicy() {}

  explicit GuidedCoarseningPolicy(const HypernodeID) {}

  GuidedCoarseningPolicy(const GuidedCoarseningPolicy&) = delete;
  GuidedCoarseningPolicy(GuidedCoarseningPolicy&&) = delete;
  GuidedCoarseningPolicy & operator= (const GuidedCoarseningPolicy &) = delete;
  GuidedCoarseningPolicy & operator= (GuidedCoarseningPolicy &&) = delete;

  template<typename Hypergraph>
  void initialize(const Hypergraph&, const Context&, utils::Timer&) { }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool acceptContraction(const Hypergraph&, const Context&, HypernodeID, HypernodeID) const {
    return true;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void accumulate(const Context& context, EdgeMetadata& sum, const EdgeMetadata value, const HyperedgeWeight edge_weight) const {
    switch (context.coarsening.rating.ge_accumulation) {
      case GuidedEdgeAccumulation::linear:
        sum += value;
        break;
      case GuidedEdgeAccumulation::quadratic:
        sum += value / static_cast<double>(edge_weight) * value;
        break;
      case GuidedEdgeAccumulation::max:
        sum = std::max(sum, value / static_cast<float>(edge_weight));
        break;
      case GuidedEdgeAccumulation::UNDEFINED:
        // ...
        break;
    }
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool acceptEdgeContraction(const Hypergraph&,
                                                                const Context& context,
                                                                const double guiding_threshold,
                                                                const HyperedgeWeight summed_rating,
                                                                const EdgeMetadata summed_md) const {
    return computeRelativeValue(context, summed_rating, summed_md) <= guiding_threshold;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HyperedgeWeight scaledRating(const Hypergraph&,
                                                                  const Context& context,
                                                                  const double guiding_threshold,
                                                                  const HyperedgeWeight summed_rating,
                                                                  const EdgeMetadata summed_md) const {
    double scale = std::max(guiding_threshold - computeRelativeValue(context, summed_rating, summed_md), 0.0) / guiding_threshold;
    switch (context.coarsening.rating.ge_scaling) {
      case GuidedEdgeScaling::none:
        return summed_rating;
      case GuidedEdgeScaling::linear:
        return std::round(10 * scale * summed_rating);
      case GuidedEdgeScaling::quadratic:
        return std::round(10 * scale * scale * summed_rating);
      case GuidedEdgeScaling::cubic:
        return std::round(10 * scale * scale * scale * summed_rating);
      case GuidedEdgeScaling::UNDEFINED:
        // ...
        break;
    }
    return summed_rating;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE double computeRelativeValue(const Context& context,
                                                                 const HyperedgeWeight summed_rating,
                                                                 const EdgeMetadata summed_md) const {
    switch (context.coarsening.rating.ge_accumulation) {
      case GuidedEdgeAccumulation::linear:
        return summed_md / static_cast<double>(summed_rating);
      case GuidedEdgeAccumulation::quadratic:
        return std::sqrt(summed_md / static_cast<double>(summed_rating));
      case GuidedEdgeAccumulation::max:
        return summed_md;
      case GuidedEdgeAccumulation::UNDEFINED:
        // ...
        break;
    }
    return 1.0;
  }

  template<typename Hypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE double weightRatioForNode(const Hypergraph&,
                                                               const HypernodeID) const {
    return 1.0;
  }
};


using DegreeSimilarityPolicies = kahypar::meta::Typelist<PreserveRebalancingNodesPolicy, GuidedCoarseningPolicy>;

}  // namespace mt_kahypar

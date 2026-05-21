/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2026 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <functional>

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"

namespace mt_kahypar {
namespace metrics {

template<typename PartitionedHypergraph>
struct MetricsTracker {
  HyperedgeWeight objective_delta;
  vec<HypernodeWeight> part_weights;
  size_t num_overloaded;
  // TODO: ignores removed d0 nodes...
  size_t num_empty;
  double imbalance;

  const PartitionedHypergraph* phg;
  const Context* context;
  const HypernodeWeight* max_part_weights;

  MetricsTracker(const PartitionedHypergraph& phg,
                 const Context& context) :
    objective_delta(0),
    part_weights(),
    num_overloaded(0),
    num_empty(0),
    imbalance(-1.0),
    phg(&phg),
    context(&context),
    max_part_weights(context.partition.max_part_weights.data()) {
      part_weights.resize(context.partition.k);
      for (size_t i = 0; i < part_weights.size(); ++i) {
        part_weights[i] = phg.partWeight(i);
      }
      initializeConstraints();
  }

  MetricsTracker(const PartitionedHypergraph& phg,
                 const Context& context,
                 const vec<HypernodeWeight>& part_weights,
                 const std::vector<HypernodeWeight>& max_part_weights) :
    objective_delta(0),
    part_weights(part_weights),
    num_overloaded(0),
    num_empty(0),
    imbalance(-1.0),
    phg(&phg),
    context(&context),
    max_part_weights(max_part_weights.data()) {
      initializeConstraints();
  }

  void initializeConstraints() {
    num_overloaded = 0;
    num_empty = 0;
    for (size_t i = 0; i < part_weights.size(); ++i) {
      if (part_weights[i] > max_part_weights[i]) {
        num_overloaded++;
      } else if (part_weights[i] == 0) {
        num_empty++;
      }
    }
  }

  void applyMove(const Move& m) {
    ASSERT(m.isValid());
    applyMove(m.from, m.to, m.node, m.gain);
  }

  void undoMove(const Move& m) {
    ASSERT(m.isValid());
    applyMove(m.to, m.from, m.node, -m.gain);
  }

  void applyMove(PartitionID from, PartitionID to, HypernodeID node, Gain gain) {
    const HypernodeWeight node_weight = phg->nodeWeight(node);
    objective_delta -= gain;
    imbalance = -1; // invalidate

    const bool from_overloaded = part_weights[from] > max_part_weights[from];
    part_weights[from] -= node_weight;
    ASSERT(part_weights[from] >= 0);
    if (from_overloaded && part_weights[from] <= max_part_weights[from]) {
      num_overloaded--;
    }
    if (part_weights[from] == 0 && node_weight > 0) {
      num_empty++;
    }
    if (part_weights[to] == 0 && node_weight > 0) {
      num_empty--;
    }
    const bool to_overloaded = part_weights[to] > max_part_weights[to];
    part_weights[to] += node_weight;
    if (!to_overloaded && part_weights[to] > max_part_weights[to]) {
      num_overloaded++;
    }

    checkConstraints();
  }

  bool isValid() const {
    return num_overloaded == 0 && (context->partition.allow_empty_blocks || num_empty == 0);
  }

  Metrics getMetrics() {
    checkConstraints();

    return Metrics{objective_delta,
      BalanceMetrics{
        computeImbalance(),
        num_overloaded > 0,
        context->partition.allow_empty_blocks && num_empty > 0}};
  }

  bool isBetter(const Metrics& metrics) {
    // the case distinctions are not necessary, but allow us to avoid
    // recomputing the imbalance in most cases
    bool is_valid = isValid();
    bool other_is_valid = metrics.imbalance.isValidPartition();
    if (is_valid && other_is_valid) {
      if (objective_delta < metrics.quality) {
        return true;
      } else if (objective_delta == metrics.quality) {
        return getMetrics().isBetter(metrics);
      } else {
        return false;
      }
    } else if (is_valid && !other_is_valid) {
      return true;
    } else if (!is_valid && other_is_valid) {
      return false;
    } else {
      return getMetrics().isBetter(metrics);
    }
  }

  // ! for parallel prefix sum
  MetricsTracker splitPrescan() const {
    MetricsTracker copy(*this);
    copy.part_weights.assign(part_weights.size(), 0);
    copy.objective_delta = 0;
    copy.imbalance = -1;
    return copy;
  }

  // ! for parallel prefix sum
  void applyMovePreScan(const Move& m) {
    ASSERT(m.isValid() && imbalance == -1);
    const HypernodeWeight node_weight = phg->nodeWeight(m.node);
    part_weights[m.from] -= node_weight;
    part_weights[m.to] += node_weight;
    objective_delta -= m.gain;
  }

  // ! for parallel prefix sum
  void addPrefix(MetricsTracker& lhs) {
    for (size_t i = 0; i < part_weights.size(); ++i) {
      part_weights[i] += lhs.part_weights[i];
    }
    objective_delta += lhs.objective_delta;
  }

 private:
  double computeImbalance() {
    if (imbalance == -1) {
      for (size_t i = 0; i < part_weights.size(); ++i) {
        const double balance_i = (part_weights[i]
                / static_cast<double>(context->partition.perfect_balance_part_weights[i]));
        imbalance = std::max(imbalance, balance_i);
      }
      imbalance -= 1;
    }
    return imbalance;
  }

  void checkConstraints() const {
    ASSERT([&]{
      MetricsTracker copy(*this);
      copy.initializeConstraints();
      if (num_overloaded != copy.num_overloaded) {
        LOG << V(num_overloaded) << V(copy.num_overloaded);
        return false;
      }
      if (num_empty != copy.num_empty) {
        LOG << V(num_empty) << V(copy.num_empty);
        return false;
      }
      return true;
    }());
  }
};

}  // namespace metrics
}  // namespace mt_kahypar

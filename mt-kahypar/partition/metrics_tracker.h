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

template<typename PartitionedHypergraph, typename Subclass>
struct MetricsTrackerBase {
  HyperedgeWeight objective_delta;
  vec<HypernodeWeight> part_weights;
  size_t num_overloaded;
  // TODO: ignores removed d0 nodes...
  size_t num_underloaded;
  double imbalance;

  const PartitionedHypergraph* phg;
  const Context* context;
  const HypernodeWeight* max_part_weights;

  MetricsTrackerBase(const PartitionedHypergraph& phg,
                     const Context& context) :
    objective_delta(0),
    part_weights(),
    num_overloaded(0),
    num_underloaded(0),
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

  MetricsTrackerBase(const PartitionedHypergraph& phg,
                     const Context& context,
                     const vec<HypernodeWeight>& part_weights,
                     const std::vector<HypernodeWeight>& max_part_weights) :
    objective_delta(0),
    part_weights(part_weights),
    num_overloaded(0),
    num_underloaded(0),
    imbalance(-1.0),
    phg(&phg),
    context(&context),
    max_part_weights(max_part_weights.data()) {
      initializeConstraints();
  }

  void initializeConstraints() {
    num_overloaded = 0;
    num_underloaded = 0;
    static_cast<Subclass*>(this)->clear();

    for (size_t i = 0; i < part_weights.size(); ++i) {
      if (isOverloaded(i)) {
        num_overloaded++;
        static_cast<Subclass*>(this)->setOverloaded(i, true);
      } else if (isUnderloaded(i)) {
        num_underloaded++;
        static_cast<Subclass*>(this)->setUnderloaded(i, true);
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

    const bool from_was_overloaded = isOverloaded(from);
    const bool from_was_underloaded = isUnderloaded(from);
    part_weights[from] -= node_weight;
    ASSERT(part_weights[from] >= 0);
    if (from_was_overloaded && !isOverloaded(from)) {
      num_overloaded--;
      static_cast<Subclass*>(this)->setOverloaded(from, false);
    }
    if (!from_was_underloaded && isUnderloaded(from)) {
      num_underloaded++;
      static_cast<Subclass*>(this)->setUnderloaded(from, true);
    }

    const bool to_was_overloaded = isOverloaded(to);
    const bool to_was_underloaded = isUnderloaded(to);
    part_weights[to] += node_weight;
    if (!to_was_overloaded && isOverloaded(to)) {
      num_overloaded++;
      static_cast<Subclass*>(this)->setOverloaded(to, true);
    }
    if (to_was_underloaded && !isUnderloaded(to)) {
      num_underloaded--;
      static_cast<Subclass*>(this)->setUnderloaded(to, false);
    }

    checkConstraints();
  }

  bool isValid() const {
    return num_overloaded == 0 && num_underloaded == 0;
  }

  Metrics getMetrics() {
    checkConstraints();

    BalanceMetrics balance{computeImbalance(), num_overloaded > 0, num_underloaded > 0};
    return Metrics{objective_delta, balance};
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
  Subclass splitPrescan() const {
    Subclass copy(static_cast<const Subclass&>(*this));
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
  void addPrefix(const MetricsTrackerBase& lhs) {
    for (size_t i = 0; i < part_weights.size(); ++i) {
      part_weights[i] += lhs.part_weights[i];
    }
    objective_delta += lhs.objective_delta;
  }

 private:
  bool isOverloaded(PartitionID block) {
    return part_weights[block] > max_part_weights[block];
  }

  bool isUnderloaded(PartitionID block) {
    ASSERT(part_weights[block] >= 0);
    return !context->partition.allow_empty_blocks && part_weights[block] == 0;
  }

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
      Subclass copy(static_cast<const Subclass&>(*this));
      copy.initializeConstraints();
      if (num_overloaded != copy.num_overloaded) {
        LOG << V(num_overloaded) << V(copy.num_overloaded);
        return false;
      }
      if (num_underloaded != copy.num_underloaded) {
        LOG << V(num_underloaded) << V(copy.num_underloaded);
        return false;
      }
      return true;
    }());
  }
};

template<typename PartitionedHypergraph>
struct MetricsTracker: public MetricsTrackerBase<PartitionedHypergraph, MetricsTracker<PartitionedHypergraph>> {
  using Base = MetricsTrackerBase<PartitionedHypergraph, MetricsTracker<PartitionedHypergraph>>;

  MetricsTracker(const PartitionedHypergraph& phg,
                 const Context& context) :
    Base(phg, context) { }

  MetricsTracker(const PartitionedHypergraph& phg,
                 const Context& context,
                 const vec<HypernodeWeight>& part_weights,
                 const std::vector<HypernodeWeight>& max_part_weights) :
    Base(phg, context, part_weights, max_part_weights) { }

  void clear() { }
  void setOverloaded(PartitionID, bool) { }
  void setUnderloaded(PartitionID, bool) { }
};

template<typename PartitionedHypergraph>
struct MetricsAndBlockTracker: public MetricsTrackerBase<PartitionedHypergraph, MetricsAndBlockTracker<PartitionedHypergraph>> {
  using Base = MetricsTrackerBase<PartitionedHypergraph, MetricsAndBlockTracker<PartitionedHypergraph>>;

  vec<PartitionID> overloaded_blocks;
  vec<PartitionID> underloaded_blocks;

  MetricsAndBlockTracker(const PartitionedHypergraph& phg,
                         const Context& context) :
    Base(phg, context) { }

  MetricsAndBlockTracker(const PartitionedHypergraph& phg,
                         const Context& context,
                         const vec<HypernodeWeight>& part_weights,
                         const std::vector<HypernodeWeight>& max_part_weights) :
    Base(phg, context, part_weights, max_part_weights) { }

  void clear() {
    overloaded_blocks.clear();
    underloaded_blocks.clear();
  }

  void setOverloaded(PartitionID block, bool value) {
    setBlockImpl(overloaded_blocks, block, value);
  }

  void setUnderloaded(PartitionID block, bool value) {
    setBlockImpl(underloaded_blocks, block, value);
  }

 private:
  void setBlockImpl(vec<PartitionID>& block_set, PartitionID block, bool value) {
    ASSERT([&]{
      bool contained = false;
      for (PartitionID p: block_set) {
        contained |= (p == block);
      }
      return contained != value;
    }());

    if (value) {
      block_set.push_back(block);
    } else {
      block_set.erase(
        std::remove_if(block_set.begin(), block_set.end(),
            [block](PartitionID p) { return p == block; }),
        block_set.end());
    }
  }
};

}  // namespace metrics
}  // namespace mt_kahypar

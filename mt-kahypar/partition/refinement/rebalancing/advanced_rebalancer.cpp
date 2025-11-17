/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/refinement/rebalancing/advanced_rebalancer.h"

#include <array>
#include <optional>
#include <random>

#include <boost/range/irange.hpp>

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

namespace impl {
  double imbalance(const HypernodeWeightArray& part_weights, const Context& context) {
    double max_balance = 0;
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      for (Dimension d = 0; d < part_weights.dimension(); ++d) {
        const double curr_balance =
                (part_weights[i].at(d) /
                  static_cast<double>(context.partition.perfect_balance_part_weights[i].at(d)));
        max_balance = std::max(max_balance, curr_balance);
      }
    }

    return max_balance - 1.0;
  }

  double imbalanceSum(const HypernodeWeightArray& part_weights, const Context& context, const vec<float>& weight_normalizer) {
    double sum = 0;
    for (PartitionID i = 0; i < context.partition.k; ++i) {
      for (Dimension d = 0; d < part_weights.dimension(); ++d) {
        HNWeightScalar diff = part_weights[i].at(d) - context.partition.max_part_weights[i].at(d);
        if (diff > 0) {
          sum += weight_normalizer[d] * diff;
        }
      }
    }
    return sum;
  }

  float transformGainFromProgress(Gain gain_, float progress, bool has_negative_progress, double negative_progress_penalty) {
    // here: positive gain means improvement
    float gain = gain_;
    if (has_negative_progress) {
      progress /= negative_progress_penalty;
    }
    if (gain > 0) {
      gain *= progress;
    } else if (gain < 0) {
      gain /= progress;
    }
    return gain;
  }

  float transformGain(Gain gain_, HNWeightConstRef wu, HNWeightAtomicCRef from_weight, HNWeightConstRef max_part_weight) {
    // here: positive gain means improvement
    if (wu.dimension() == 1) {
      float gain = gain_;
      if (gain > 0) {
        gain *= wu.at(0);
      } else if (gain < 0) {
        gain /= wu.at(0);
      }
      return gain;
    } else {
      float relevant_weight_fraction = 0;
      for (Dimension d = 0; d < wu.dimension(); ++d) {
        if (from_weight.at(d) > max_part_weight.at(d)) {
          relevant_weight_fraction += wu.at(d) / static_cast<float>(max_part_weight.at(d));
        }
      }
      return transformGainFromProgress(gain_, relevant_weight_fraction, false, 1.0);
    }
  }

  std::pair<float, bool> computeBalanceProgress(HNWeightConstRef wu, HNWeightConstRef from_weight, HNWeightConstRef max_part_weight_from,
                                                HNWeightAtomicCRef to_weight, HNWeightConstRef max_part_weight_to, const vec<float>& weight_normalizer) {
    ASSERT(wu.dimension() > 1);
    float relative_progress = 0;
    const auto new_to_weight = wu + to_weight;
    bool has_negative_progress = false;
    for (Dimension d = 0; d < wu.dimension(); ++d) {
      if (from_weight.at(d) > max_part_weight_from.at(d)) {
        relative_progress += weight_normalizer[d] * wu.at(d);
      }
      const HNWeightScalar overweight = std::min(new_to_weight.at(d) - max_part_weight_to.at(d), wu.at(d));
      if (overweight > 0) {
        relative_progress -= weight_normalizer[d] * overweight;
        has_negative_progress = true;
      }
    }
    return {relative_progress, has_negative_progress};
  }

  template<typename PartitionedHypergraph, typename GainCache, typename Range>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  std::pair<PartitionID, float> computeBestTargetBlockGeneric(
          const PartitionedHypergraph& phg, const Context& context, const GainCache& gain_cache,
          HypernodeID u, PartitionID from, Range range,
          const HypernodeWeightArray& reduced_part_weights, bool skip_non_adjacent,
          AllocatedHNWeight& best_to_weight, AllocatedHNWeight& tmp_hn_weight,
          const vec<float>& weight_normalizer) {
    const HNWeightConstRef wu = phg.nodeWeight(u);
    const HNWeightAtomicCRef from_weight = phg.partWeight(from);
    const bool any_progress = context.refinement.rebalancing.allow_any_progress;

    if (!any_progress) {
      // the normal case
      PartitionID to = kInvalidPartition;
      HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
      best_to_weight = from_weight - wu;
      for (PartitionID i : range) {
        if (i != from && i != kInvalidPartition) {
          tmp_hn_weight = phg.partWeight(i);
          HyperedgeWeight benefit;
          if (gain_cache.blockIsAdjacent(u, i)) {
            benefit = gain_cache.benefitTerm(u, i);
          } else if (skip_non_adjacent && (to != kInvalidPartition || !(tmp_hn_weight + wu <= reduced_part_weights[i])) ) {
            // skip expensive gain recomputation
            continue;
          } else {
            benefit = gain_cache.recomputeBenefitTerm(phg, u, i);
          }
          // TODO: any better tie breaking option?
          if ((benefit > to_benefit || (benefit == to_benefit && tmp_hn_weight < best_to_weight)) &&
              tmp_hn_weight + wu <= reduced_part_weights[i]) {
            to_benefit = benefit;
            to = i;
            best_to_weight = tmp_hn_weight;
          }
        }
      }

      if (to != kInvalidPartition) {
        Gain gain = to_benefit - gain_cache.penaltyTerm(u, phg.partID(u));
        return std::make_pair(to, transformGain(gain, wu, phg.partWeight(from), reduced_part_weights[from]));
      }
    } else {
      // here, we need to compare the gains with included weight instead fo simply the move gain
      best_to_weight = phg.partWeight(from);  // abuse this for caching the current part weight
      const HNWeightConstRef max_from_weight = reduced_part_weights[from];

      PartitionID to = kInvalidPartition;
      HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
      float to_gain = std::numeric_limits<float>::min();
      for (PartitionID i : range) {
        if (i != from && i != kInvalidPartition) {
          bool is_adjacent;
          if (skip_non_adjacent) {
            is_adjacent = gain_cache.blockIsAdjacent(u, i);
            if (!is_adjacent && to != kInvalidPartition) continue;  // skip expensive gain recomputation
          }

          auto [progress, negative_progress] = computeBalanceProgress(wu, best_to_weight, max_from_weight,
                                                                      phg.partWeight(i), reduced_part_weights[i], weight_normalizer);
          if (progress > 0) {
            if (!skip_non_adjacent) {
              is_adjacent = gain_cache.blockIsAdjacent(u, i);
            }
            HyperedgeWeight benefit = is_adjacent ? gain_cache.benefitTerm(u, i) : gain_cache.recomputeBenefitTerm(phg, u, i);
            float gain = transformGainFromProgress(benefit, progress, negative_progress,
                                                   context.refinement.rebalancing.negative_progress_penalty);
            if (gain > to_gain) {
              to_gain = gain;
              to_benefit = benefit;
              to = i;
            }
          }
        }
      }

      if (to != kInvalidPartition) {
        Gain gain = to_benefit - gain_cache.penaltyTerm(u, phg.partID(u));
        best_to_weight = phg.partWeight(from);
        auto [progress, negative_progress] = computeBalanceProgress(wu, best_to_weight, max_from_weight,
                                                                    phg.partWeight(to), reduced_part_weights[to], weight_normalizer);
        return std::make_pair(to, transformGainFromProgress(gain, progress, negative_progress,
                                                            context.refinement.rebalancing.negative_progress_penalty));
      }
    }
    return std::make_pair(kInvalidPartition, std::numeric_limits<float>::min());
  }

  template<typename PartitionedHypergraph, typename GainCache>
  std::pair<PartitionID, float> computeBestTargetBlock(
          const PartitionedHypergraph& phg, const Context& context, const GainCache& gain_cache,
          HypernodeID u, PartitionID from, const HypernodeWeightArray& reduced_part_weights,
          AllocatedHNWeight& best_to_weight, AllocatedHNWeight& tmp_hn_weight, const vec<float>& weight_normalizer) {
    return computeBestTargetBlockGeneric(phg, context, gain_cache, u, from, boost::irange<PartitionID>(0, context.partition.k),
                                         reduced_part_weights, true, best_to_weight, tmp_hn_weight, weight_normalizer);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  std::pair<PartitionID, float> bestOfThree(
          const PartitionedHypergraph& phg, const Context& context, const GainCache& gain_cache,
          HypernodeID u, PartitionID from, std::array<PartitionID, 3> parts, const HypernodeWeightArray& reduced_part_weights,
          AllocatedHNWeight& best_to_weight, AllocatedHNWeight& tmp_hn_weight, const vec<float>& weight_normalizer) {
    auto [to, gain] = computeBestTargetBlockGeneric(phg, context, gain_cache, u, from, IteratorRange(parts.cbegin(), parts.cend()),
                                                    reduced_part_weights, false, best_to_weight, tmp_hn_weight, weight_normalizer);

    if (to != kInvalidPartition) {
      return std::make_pair(to, gain);
    } else {
      // edge case: if u does not fit in any of the three considered blocks we need to check all blocks
      return computeBestTargetBlock(phg, context, gain_cache, u, from, reduced_part_weights, best_to_weight, tmp_hn_weight, weight_normalizer);
    }
  }

  struct AccessToken {
    AccessToken(int seed, size_t num_pqs) : dist(0, num_pqs - 1) { rng.seed(seed); }
    size_t getRandomPQ() { return dist(rng); }

    std::array<size_t, 2> getTwoRandomPQs() {
      std::array<size_t, 2> result({getRandomPQ(), getRandomPQ()});
      while (result[0] == result[1]) { result[1] = getRandomPQ(); }
      return result;
    }

    std::mt19937 rng;
    std::uniform_int_distribution<size_t> dist;
  };


  template<typename PartitionedHypergraph, typename GainCache>
  struct NextMoveFinder {
    Move next_move;

    PartitionedHypergraph& _phg;
    GainCache& _gain_cache;
    const Context& _context;
    const HypernodeWeightArray& _reduced_part_weights;
    const vec<float>& _weight_normalizer;

    vec<rebalancer::GuardedPQ>& _pqs;
    ds::Array<PartitionID>& _target_part;
    ds::Array<rebalancer::NodeState>& _node_state;
    AllocatedHNWeight& _best_to_weight;
    AllocatedHNWeight& _tmp_hn_weight;
    AccessToken _token;

    NextMoveFinder(int seed, const Context& context, const HypernodeWeightArray& reduced_part_weights, const vec<float>& weight_normalizer,
                   PartitionedHypergraph& phg, GainCache& gain_cache, vec<rebalancer::GuardedPQ>& pqs,
                   ds::Array<PartitionID>& target_part, ds::Array<rebalancer::NodeState>& node_state,
                   AllocatedHNWeight& best_to_weight, AllocatedHNWeight& tmp_hn_weight) :
                   _phg(phg), _gain_cache(gain_cache), _context(context), _reduced_part_weights(reduced_part_weights), _weight_normalizer(weight_normalizer),
                   _pqs(pqs), _target_part(target_part), _node_state(node_state),
                   _best_to_weight(best_to_weight), _tmp_hn_weight(tmp_hn_weight), _token(seed, pqs.size()) { }


    void recomputeTopGainMove(HypernodeID v, const Move& move /* of the neighbor */) {
      float gain = 0;
      PartitionID newTarget = kInvalidPartition;
      const PartitionID designatedTargetV = _target_part[v];
      if (_context.partition.k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
        std::tie(newTarget, gain) = computeBestTargetBlock(_phg, _context, _gain_cache, v, _phg.partID(v),
                                                           _reduced_part_weights, _best_to_weight, _tmp_hn_weight, _weight_normalizer);
      } else {
        std::tie(newTarget, gain) = bestOfThree(_phg, _context, _gain_cache,
                                                v, _phg.partID(v), {designatedTargetV, move.from, move.to},
                                                _reduced_part_weights, _best_to_weight, _tmp_hn_weight, _weight_normalizer);
      }
      _target_part[v] = newTarget;
    }

    bool checkCandidate(HypernodeID u, float& gain_in_pq) {
      if (!_node_state[u].tryLock()) return false;
      auto [to, true_gain] = computeBestTargetBlock(_phg, _context, _gain_cache, u, _phg.partID(u),
                                                    _reduced_part_weights, _best_to_weight, _tmp_hn_weight, _weight_normalizer);
      if (to != kInvalidPartition && true_gain >= gain_in_pq) {
        next_move.node = u;
        next_move.to = to;
        next_move.from = _phg.partID(u);
        next_move.gain = true_gain;
        return true;
      } else {
        _target_part[u] = to;
        gain_in_pq = true_gain;
        _node_state[u].unlock();
        return false;
      }
    }

    bool lockedModifyPQ(size_t best_id, bool parallel) {
      auto& gpq = _pqs[best_id];
      auto& pq = gpq.pq;

      HypernodeID node = pq.top();
      float gain_in_pq = pq.topKey();
      const bool success = checkCandidate(node, gain_in_pq);

      if (success) {
        pq.deleteTop();
        gpq.top_key = pq.empty() ? std::numeric_limits<float>::lowest() : pq.topKey();
      } else {
        // gain was updated by success_func in this case
        if (_target_part[node] != kInvalidPartition) {
          pq.adjustKey(node, gain_in_pq);
          gpq.top_key = pq.topKey();
        } else {
          pq.deleteTop();
          gpq.top_key = pq.empty() ? std::numeric_limits<float>::lowest() : pq.topKey();
        }
      }
      if (parallel) {
        gpq.lock.unlock();
      }
      return success;
    }

    bool tryPop() {
      static constexpr size_t NUM_TRIES = 32;
      for (size_t i = 0; i < NUM_TRIES; ++i) {
        auto two = _token.getTwoRandomPQs();
        auto& first = _pqs[two[0]];
        auto& second = _pqs[two[1]];
        if (first.pq.empty() && second.pq.empty()) continue;
        size_t best_id = two[0];
        if (first.pq.empty() || first.top_key < second.top_key) best_id = two[1];
        if (!_pqs[best_id].lock.tryLock()) continue;
        // could also check for top key. would want to distinguish tries that failed due to high contention
        // vs approaching the end
        if (_pqs[best_id].pq.empty()) {
          _pqs[best_id].lock.unlock();
          continue;
        }
        if (lockedModifyPQ(best_id, true)) return true;
        // if you got a PQ but it fails because the node's gain was wrong or the node couldn't be locked
        // (success_func failed) then we still want to use the standard method
        i = 0;
      }

      while (true) {
        float best_key = std::numeric_limits<float>::lowest();
        int best_id = -1;
        for (size_t i = 0; i < _pqs.size(); ++i) {
          if (!_pqs[i].pq.empty() && _pqs[i].top_key > best_key) {
            best_key = _pqs[i].top_key;
            best_id = i;
          }
        }
        if (best_id == -1) return false;
        if (!_pqs[best_id].lock.tryLock()) continue;
        if (_pqs[best_id].pq.empty()) {
          _pqs[best_id].lock.unlock();
          continue;
        }
        if (lockedModifyPQ(best_id, true)) return true;
      }
    }

    bool tryPopSequential() {
      while (true) {
        float best_key = std::numeric_limits<float>::min();
        int best_id = -1;
        for (size_t i = 0; i < _pqs.size(); ++i) {
          if (!_pqs[i].pq.empty() && _pqs[i].top_key > best_key) {
            best_key = _pqs[i].top_key;
            best_id = i;
          }
        }
        if (best_id == -1) return false;
        ASSERT(!_pqs[best_id].pq.empty());
        if (lockedModifyPQ(best_id, false)) return true;
      }
    }

    bool findNextMove(bool parallel) {
      if (parallel) {
        return tryPop();
      } else {
        return tryPopSequential();
      }
    }
  };

  void deactivateOverloadedBlock(uint8_t* is_overloaded, size_t* num_overloaded_blocks) {
    if (*is_overloaded) {
      uint8_t expected = 1;
      if (__atomic_compare_exchange_n(is_overloaded, &expected, 0, false, __ATOMIC_ACQUIRE, __ATOMIC_RELAXED)) {
        __atomic_fetch_sub(num_overloaded_blocks, 1, __ATOMIC_RELAXED);
      }
    }
  }

} // namespace impl


  static constexpr MoveID kInvalidMove = std::numeric_limits<MoveID>::max();

  template <typename GraphAndGainTypes>
  void AdvancedRebalancer<GraphAndGainTypes>::insertNodesInOverloadedBlocks(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                            const HypernodeWeightArray& reduced_part_weights) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    // init PQs if not done before
    const size_t num_pqs = 2 * _context.shared_memory.num_threads;
    if (_pqs.size() != num_pqs) {
      _pqs.assign(num_pqs, rebalancer::GuardedPQ(_pq_handles.data(), _node_state.size()));
    }
    for (auto& gpq : _pqs) {
      gpq.reset();
    }

    // data structures to draw random PQs
    std::atomic<int> seed { 555 };
    tbb::enumerable_thread_specific<impl::AccessToken> ets_tokens([&]() {
      return impl::AccessToken(seed.fetch_add(1, std::memory_order_relaxed), num_pqs);
    });

    // insert nodes into PQs
    phg.doParallelForAllNodes([&](HypernodeID u) {
      const PartitionID b = phg.partID(u);
      if (!_is_overloaded[b] || phg.isFixed(u)) return;

      auto [target, gain] = impl::computeBestTargetBlock(phg, _context, _gain_cache, u, phg.partID(u), reduced_part_weights,
                                                         _best_target_block_weight.local(), _tmp_hn_weight.local(), _weight_normalizer);
      ASSERT(target == kInvalidPartition || _context.refinement.rebalancing.allow_any_progress ||
             gain == impl::transformGain(_gain_cache.recomputeBenefitTerm(phg, u, target) - _gain_cache.recomputePenaltyTerm(phg, u),
                                         phg.nodeWeight(u), phg.partWeight(phg.partID(u)), reduced_part_weights[phg.partID(u)]),
             "Gain cache is in invalid state!");
      if (target == kInvalidPartition) return;

      _node_state[u].markAsMovable();
      _target_part[u] = target;

      auto& token = ets_tokens.local();
      int my_pq_id = -1;
      while (true) {
        my_pq_id = token.getRandomPQ();
        if (_pqs[my_pq_id].lock.tryLock()) {
          break;
        }
      }
      _pqs[my_pq_id].pq.insert(u, gain);
      _pqs[my_pq_id].lock.unlock();
      _pq_id[u] = my_pq_id;
    });


    for (rebalancer::GuardedPQ& gpq : _pqs) {
      if (!gpq.pq.empty()) {
        gpq.top_key = gpq.pq.topKey();
      }
    }
  }

  template <typename GraphAndGainTypes>
  int64_t AdvancedRebalancer<GraphAndGainTypes>::findMoves(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                           const HypernodeWeightArray& reduced_part_weights,
                                                           size_t& global_move_id,
                                                           bool parallel) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    const bool any_progress = _context.refinement.rebalancing.allow_any_progress;
    int64_t attributed_gain = 0;
    size_t num_overloaded_blocks = _overloaded_blocks.size();

    AllocatedHNWeight inf_weight;
    if (any_progress) {
      inf_weight = weight::broadcast(std::numeric_limits<HNWeightScalar>::max(), phg.dimension());
    }

    auto task = [&](size_t task_id) {
      vec<HyperedgeID> edges_with_gain_changes;
      Gain local_attributed_gain = 0;
      vec<vec<HypernodeID>> nodes_to_update(_pqs.size());
      vec<int> pqs_to_update;

      const int seed = phg.initialNumNodes() + task_id;

      impl::NextMoveFinder<PartitionedHypergraph, GainCache> next_move_finder(
        seed, _context, reduced_part_weights, _weight_normalizer, phg, _gain_cache, _pqs, _target_part, _node_state,
        _best_target_block_weight.local(), _tmp_hn_weight.local());

      while (num_overloaded_blocks > 0 && next_move_finder.findNextMove(parallel)) {
        const Move& m = next_move_finder.next_move;
        ASSERT(m.to != kInvalidPartition);
        const PartitionID from = phg.partID(m.node);
        _node_state[m.node].markAsMovedAndUnlock();

        if (phg.partWeight(from) <= _context.partition.max_part_weights[from]) {
          impl::deactivateOverloadedBlock(&_is_overloaded[from], &num_overloaded_blocks);
          continue;
        }

        edges_with_gain_changes.clear();
        size_t move_id = 0;
        bool moved = phg.changeNodePart(
                      _gain_cache, m.node, m.from, m.to,
                      any_progress ? inf_weight : _context.partition.max_part_weights[m.to],
                      [&] { move_id = __atomic_fetch_add(&global_move_id, 1, __ATOMIC_RELAXED); },
                      [&](const SynchronizedEdgeUpdate& sync_update) {
                        local_attributed_gain += AttributedGains::gain(sync_update);
                        if (!PartitionedHypergraph::is_graph && GainCache::triggersDeltaGainUpdate(sync_update)) {
                          edges_with_gain_changes.push_back(sync_update.he);
                        }
                      }
                    );



        if (!moved) continue;

        auto update_neighbor = [&](HypernodeID v) {
          if (v != m.node && _node_state[v].tryLock()) {
            int my_pq_id = _pq_id[v];
            ASSERT(my_pq_id != -1);
            if (nodes_to_update[my_pq_id].empty()) {
              pqs_to_update.push_back(my_pq_id);
            }
            nodes_to_update[my_pq_id].push_back(v);
            next_move_finder.recomputeTopGainMove(v, m);
          }
        };

        // update neighbors
        if constexpr (PartitionedHypergraph::is_graph) {
          for (const auto e : phg.incidentEdges(m.node)) {
            HypernodeID v = phg.edgeTarget(e);
            update_neighbor(v);
          }
        } else {
          for (HyperedgeID e : edges_with_gain_changes) {
            if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
              for (HypernodeID v : phg.pins(e)) {
                update_neighbor(v);
              }
            }
          }
        }

        while (!pqs_to_update.empty()) {
          for (size_t i = 0; i < pqs_to_update.size(); ++i) {
            int my_pq_id = pqs_to_update[i];
            auto& gpq = _pqs[my_pq_id];
            auto& pq = gpq.pq;
            if (gpq.lock.tryLock()) {
              for (HypernodeID v : nodes_to_update[my_pq_id]) {
                if (pq.contains(v)) {
                  if (_target_part[v] != kInvalidPartition) {
                    Gain new_gain_int;
                    const PartitionID from = phg.partID(v);
                    if (_gain_cache.blockIsAdjacent(v, _target_part[v])) {
                      new_gain_int = _gain_cache.gain(v, from, _target_part[v]);
                    } else {
                      new_gain_int = _gain_cache.recomputeBenefitTerm(phg, v, _target_part[v]) - _gain_cache.penaltyTerm(v, from);
                    }
                    float new_gain = impl::transformGain(new_gain_int, phg.nodeWeight(v), phg.partWeight(from), reduced_part_weights[from]);
                    pq.adjustKey(v, new_gain);
                  } else {
                    pq.remove(v);
                  }
                }
                _node_state[v].unlock();
              }

              gpq.lock.unlock();
              pqs_to_update[i] = pqs_to_update.back();
              pqs_to_update.pop_back();
              nodes_to_update[my_pq_id].clear();
            }
          }
        }

        ASSERT(m.isValid());
        _moves[move_id] = m;
        if (_context.refinement.rebalancing.allow_multiple_moves && _move_id_of_node[m.node] == kInvalidMove) {
          _move_id_of_node[m.node] = move_id;
        }
      }
      __atomic_fetch_add(&attributed_gain, local_attributed_gain, __ATOMIC_RELAXED);
    };

    if (parallel) {
      tbb::task_group tg;
      for (size_t i = 0; i < _context.shared_memory.num_threads; ++i) { tg.run(std::bind(task, i)); }
      tg.wait();
    } else {
      task(0);
    }

    return attributed_gain;
  }

  template <typename GraphAndGainTypes>
  int64_t AdvancedRebalancer<GraphAndGainTypes>::applyRollback(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                               const size_t old_move_id,
                                                               size_t& global_move_id) {
    // TODO: this implementation is not as performant as it could be
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    HypernodeWeightArray part_weights = phg.partWeights().copy();
    auto compute_imbalance = [&]() {
      return _context.refinement.rebalancing.l1_rollback ?
        impl::imbalanceSum(part_weights, _context, _weight_normalizer) : impl::imbalance(part_weights, _context);
    };

    double best_imbalance = compute_imbalance();
    size_t best_move_id = global_move_id;
    for (size_t i = global_move_id; i > old_move_id; --i) {
      const size_t move_id = i - 1;
      const Move& m = _moves[move_id];
      ASSERT(m.isValid());
      HNWeightConstRef weight = phg.nodeWeight(m.node);
      // revert the move
      part_weights[m.to] -= weight;
      part_weights[m.from] += weight;
      double curr_imbalance = compute_imbalance();
      // TODO: quality tie break?
      // note: < seems better than <= since it has a better chance of allowing progress
      if (curr_imbalance < best_imbalance) {
        best_imbalance = curr_imbalance;
        best_move_id = move_id;
      }
    }

    int64_t attributed_gain = 0;
    if (best_move_id < global_move_id) {
      DBG << "Rolling back" << (global_move_id - best_move_id) << "moves";
    }
    for (size_t i = global_move_id; i > best_move_id; --i) {
      const size_t move_id = i - 1;
      const Move& m = _moves[move_id];
      bool success = phg.changeNodePart(_gain_cache, m.node, m.to, m.from,
        [&](const SynchronizedEdgeUpdate& sync_update) {
          attributed_gain += AttributedGains::gain(sync_update);
        });
      ASSERT(success);
      if (_context.refinement.rebalancing.allow_multiple_moves && _move_id_of_node[m.node] == move_id) {
        _move_id_of_node[m.node] = kInvalidMove;
      }
    }
    global_move_id = best_move_id;
    return attributed_gain;
  }

  template <typename GraphAndGainTypes>
  std::pair<int64_t, size_t> AdvancedRebalancer<GraphAndGainTypes>::runGreedyRebalancingRound(
      mt_kahypar_partitioned_hypergraph_t& hypergraph,
      const HypernodeWeightArray& reduced_part_weights,
      size_t& global_move_id,
      bool parallel) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    _overloaded_blocks.clear();
    _is_overloaded.assign(_context.partition.k, false);
    for (PartitionID k = 0; k < _context.partition.k; ++k) {
      if ( !(phg.partWeight(k) <= _context.partition.max_part_weights[k]) ) {
        _overloaded_blocks.push_back(k);
        _is_overloaded[k] = 1;
      }
    }

    insertNodesInOverloadedBlocks(hypergraph, reduced_part_weights);

    const size_t old_id = global_move_id;
    int64_t attributed_gain = findMoves(hypergraph, reduced_part_weights, global_move_id, parallel);

    if (_context.refinement.rebalancing.use_rollback) {
      attributed_gain += applyRollback(hypergraph, old_id, global_move_id);
    }

    if constexpr (GainCache::invalidates_entries) {
      tbb::parallel_for(old_id, global_move_id, [&](const size_t i) {
        _gain_cache.recomputeInvalidTerms(phg, _moves[i].node);
      });
    }

    phg.doParallelForAllNodes([&](HypernodeID u) {
      _node_state[u].reset();
    });

    for (auto& gpq : _pqs) {
      gpq.pq.clear();
    }

    size_t num_overloaded_blocks = 0;
    for (PartitionID b = 0; b < _context.partition.k; ++b) {
      if ( !(phg.partWeight(b) <= _context.partition.max_part_weights[b]) ) {
        num_overloaded_blocks++;
      }
    }
    DBG << "Rebalancing round: moved" << (global_move_id - old_id) << "nodes; gain =" << attributed_gain << " new imbalance =" << metrics::imbalance(phg, _context);

    return {attributed_gain, num_overloaded_blocks};
  }

  template <typename GraphAndGainTypes>
  std::tuple<int64_t, size_t, size_t> AdvancedRebalancer<GraphAndGainTypes>::runGreedyAlgorithm(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                                                size_t& global_move_id) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    int64_t attributed_gain = 0;
    size_t num_overloaded_blocks = 0;
    double old_overweight = impl::imbalanceSum(phg.partWeights(), _context, _weight_normalizer);
    double new_overweight = old_overweight;
    size_t num_moves_first_round = 0;
    bool moved_nodes = false;
    do {
      old_overweight = new_overweight;
      const size_t old_id = global_move_id;
      auto [attr_gain, n_overloaded] =
        runGreedyRebalancingRound(hypergraph, _context.partition.max_part_weights, global_move_id, true);
      attributed_gain += attr_gain;
      num_overloaded_blocks = n_overloaded;
      new_overweight = impl::imbalanceSum(phg.partWeights(), _context, _weight_normalizer);
      moved_nodes = global_move_id > old_id;
      if (num_moves_first_round == 0) {
        num_moves_first_round = global_move_id;
      }
      ASSERT((num_overloaded_blocks == 0) == (new_overweight == 0));
    } while (_context.refinement.rebalancing.allow_multiple_moves
             && num_overloaded_blocks > 0
             && new_overweight < old_overweight
             && global_move_id < phg.initialNumNodes());
    DBG << V(old_overweight) << V(new_overweight) << V(global_move_id) << V(moved_nodes);

    if (_context.refinement.rebalancing.finalize_sequential && !moved_nodes) {
      auto [attr_gain, n_overloaded] =
        runGreedyRebalancingRound(hypergraph, _context.partition.max_part_weights, global_move_id, false);
      attributed_gain += attr_gain;
      num_overloaded_blocks = n_overloaded;
      old_overweight = new_overweight;
      new_overweight = impl::imbalanceSum(phg.partWeights(), _context, _weight_normalizer);
      DBG << "sequential round:" << V(old_overweight) << V(new_overweight) << V(global_move_id);
    }
    return {attributed_gain, num_overloaded_blocks, num_moves_first_round};
  }

  template <typename GraphAndGainTypes>
  bool AdvancedRebalancer<GraphAndGainTypes>::refineInternalParallel(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                  vec<vec<Move>>* moves_by_part,
                                                                  vec<Move>* moves_linear,
                                                                  Metrics& best_metric) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(_gain_cache));
    DBG << "Rebalancing: initial imbalance =" << best_metric.imbalance << " initial cut =" << best_metric.quality;

    if (_context.refinement.rebalancing.allow_multiple_moves) {
      _move_id_of_node.assign(phg.initialNumNodes(), kInvalidMove);
    }
    _weight_normalizer.resize(phg.dimension(), 0);
    for (Dimension d = 0; d < phg.dimension(); ++d) {
      _weight_normalizer[d] = 1 / static_cast<float>(phg.totalWeight().at(d));
    }

    size_t global_move_id = 0;
    auto [attributed_gain, num_overloaded_blocks, num_moves_first_round] = runGreedyAlgorithm(hypergraph, global_move_id);

    if (_context.refinement.rebalancing.allow_multiple_moves
        && (moves_by_part != nullptr || moves_linear != nullptr)) {
      // deduplicate moves: callers are allowed to assume that each move is is unique in the move sequence
      for (size_t i = num_moves_first_round; i < global_move_id; ++i) {
        Move& r_move = _moves[i];
        ASSERT(r_move.isValid() && _move_id_of_node[r_move.node] != kInvalidMove);
        if (_move_id_of_node[r_move.node] < i) {
          // "merge" the moves
          Move& first_move = _moves[_move_id_of_node[r_move.node]];
          ASSERT(first_move.isValid() && r_move.node == first_move.node);
          first_move.to = r_move.to;
          r_move.invalidate();
        }
      }
    }

    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation(_gain_cache));

    if (moves_by_part != nullptr) {
      moves_by_part->resize(_context.partition.k);
      for (auto& direction : *moves_by_part) direction.clear();
      for (size_t i = 0; i < global_move_id; ++i) {
        const Move& m = _moves[i];
        if (m.isValid() && m.from != m.to) {
          (*moves_by_part)[m.from].push_back(m);
        }
      }
    } else if (moves_linear != nullptr) {
      moves_linear->clear();
      moves_linear->reserve(global_move_id);
      for (size_t i = 0; i < global_move_id; ++i) {
        const Move& m = _moves[i];
        if (m.isValid() && m.from != m.to) {
          moves_linear->push_back(m);
        }
      }
    }

    best_metric.quality += attributed_gain;
    best_metric.imbalance = metrics::imbalance(phg, _context);
    if (num_overloaded_blocks > 0) {
      DBG << RED << "Rebalancing:   final imbalance =" << best_metric.imbalance << " final cut =" << best_metric.quality << END;
    } else {
      DBG << GREEN << "Rebalancing:   final imbalance =" << best_metric.imbalance << " final cut =" << best_metric.quality << END;
    }

    return num_overloaded_blocks == 0;
  }


template <typename GraphAndGainTypes>
AdvancedRebalancer<GraphAndGainTypes>::AdvancedRebalancer(
        HypernodeID num_nodes, const Context& context, GainCache& gain_cache) :
        _context(context),
        _gain_cache(gain_cache),
        _current_k(_context.partition.k),
        _gain(context),
        _moves(2 * num_nodes),
        _move_id_of_node(num_nodes),
        _target_part(num_nodes, kInvalidPartition),
        _pq_handles(num_nodes, invalid_position),
        _pq_id(num_nodes, -1),
        _node_state(num_nodes) { }

template <typename GraphAndGainTypes>
AdvancedRebalancer<GraphAndGainTypes>::AdvancedRebalancer(
        HypernodeID num_nodes, const Context& context, gain_cache_t gain_cache) :
        AdvancedRebalancer(num_nodes, context, GainCachePtr::cast<GainCache>(gain_cache)) { }


template <typename GraphAndGainTypes>
bool AdvancedRebalancer<GraphAndGainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                const vec<HypernodeID>& , Metrics& best_metrics, double) {
  return refineInternalParallel(hypergraph, nullptr, nullptr, best_metrics);
}

template <typename GraphAndGainTypes>
void AdvancedRebalancer<GraphAndGainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
  auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

  if (!_gain_cache.isInitialized()) {
    _gain_cache.initializeGainCache(phg);
  }
}

template <typename GraphAndGainTypes>
bool AdvancedRebalancer<GraphAndGainTypes>::refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                  const vec<HypernodeID>& ,
                                                                  vec<vec<Move>>& moves_by_part,
                                                                  Metrics& best_metrics,
                                                                  const double) {
  return refineInternalParallel(hypergraph, &moves_by_part, nullptr, best_metrics);
}

template <typename GraphAndGainTypes>
bool AdvancedRebalancer<GraphAndGainTypes>::refineAndOutputMovesLinearImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                        const vec<HypernodeID>& ,
                                                                        vec<Move>& moves,
                                                                        Metrics& best_metrics,
                                                                        const double) {
  return refineInternalParallel(hypergraph, nullptr, &moves, best_metrics);
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
  #define ADVANCED_REBALANCER(X) AdvancedRebalancer<X>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_VALID_TRAITS(ADVANCED_REBALANCER)

}   // namespace mt_kahypar

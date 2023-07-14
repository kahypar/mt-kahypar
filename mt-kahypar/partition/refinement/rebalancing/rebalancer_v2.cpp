#include <mt-kahypar/datastructures/priority_queue.h>
#include <optional>
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_v2.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/partition/context.h"

#include "pcg_random.hpp"

namespace mt_kahypar {

namespace impl {

  float transformGain(Gain gain, HypernodeWeight wu) {
    if (gain > 0) {
      gain *= wu;
    } else if (gain < 0) {
      gain /= wu;
    }
    return gain;
  }

  template<typename PartitionedHypergraph, typename GainCache>
  std::pair<PartitionID, float> computeBestTargetBlock(
          const PartitionedHypergraph& phg, const Context& context, const GainCache& gain_cache,
          HypernodeID u, PartitionID from) {
    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i = 0; i < phg.k(); ++i) {
      if (i != from) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight benefit = gain_cache.benefitTerm(u, i);
        if ((benefit > to_benefit || (benefit == to_benefit && to_weight < best_to_weight)) &&
        //  TODO swap this for relaxed max part weights?
            to_weight + wu <= context.partition.max_part_weights[i]) {
          to_benefit = benefit;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }

    Gain gain = std::numeric_limits<Gain>::min();
    if (to != kInvalidPartition) {
      gain = to_benefit - gain_cache.penaltyTerm(u, phg.partID(u));
    }
    return std::make_pair(to, transformGain(gain, wu));
  }

  template<typename PartitionedHypergraph, typename GainCache>
  std::pair<PartitionID, float> bestOfThree(
          const PartitionedHypergraph& phg, const Context& context, const GainCache& gain_cache,
          HypernodeID u, PartitionID from, std::array<PartitionID, 3> parts) {
    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i : parts) {
      if (i != from && i != kInvalidPartition) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight benefit = gain_cache.benefitTerm(u, i);
        if ((benefit > to_benefit || (benefit == to_benefit && to_weight < best_to_weight)) &&
            to_weight + wu <= context.partition.max_part_weights[i]) {
          to_benefit = benefit;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }

    Gain gain = std::numeric_limits<Gain>::min();
    if (to != kInvalidPartition) {
      gain = to_benefit - gain_cache.penaltyTerm(u, phg.partID(u));
    }
    return std::make_pair(to, transformGain(gain, wu));
  }

  struct AccessToken {
    AccessToken(int seed, size_t num_pqs) : dist(0, num_pqs - 1) { rng.seed(seed); }
    size_t getRandomPQ() { return dist(rng); }

    std::array<size_t, 2> getTwoRandomPQs() {
      std::array<size_t, 2> result({getRandomPQ(), getRandomPQ()});
      while (result[0] != result[1]) { result[1] = getRandomPQ(); }
      return result;
    }

    pcg32 rng;
    std::uniform_int_distribution<size_t> dist;
  };


  template<typename PartitionedHypergraph, typename GainCache>
  struct NextMoveFinder {

    NextMoveFinder(int seed, const Context& context, PartitionedHypergraph& phg, GainCache& gain_cache,
                   vec<rebalancer::GuardedPQ>& pqs, vec<PartitionID>& target_part, vec<rebalancer::NodeState>& node_state) :
                   _phg(phg), _gain_cache(gain_cache), _context(context),
                   _pqs(pqs), _target_part(target_part), _node_state(node_state), _token(seed, pqs.size()) { }


    void recomputeTopGainMove(HypernodeID v, const Move& move /* of the neighbor */) {
      float gain = 0;
      PartitionID newTarget = kInvalidPartition;
      const PartitionID designatedTargetV = _target_part[v];
      if (_context.partition.k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
        std::tie(newTarget, gain) = computeBestTargetBlock(_phg, _context, _gain_cache, v, _phg.partID(v));
      } else {
        std::tie(newTarget, gain) = bestOfThree(_phg, _context, _gain_cache,
                                                v, _phg.partID(v), {designatedTargetV, move.from, move.to});
      }
      _target_part[v] = newTarget;
    }

    bool checkCandidate(HypernodeID u, float& gain_in_pq) {
      if (!_node_state[u].tryLock()) return false;
      auto [to, true_gain] = computeBestTargetBlock(_phg, _context, _gain_cache, u, _phg.partID(u));
      if (true_gain >= gain_in_pq) {
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

    bool lockedModifyPQ(size_t best_id) {
      auto& gpq = _pqs[best_id];
      auto& pq = gpq.pq;

      HypernodeID node = pq.top();
      float gain_in_pq = pq.topKey();
      const bool success = checkCandidate(node, gain_in_pq);

      if (success) {
        pq.deleteTop();
        gpq.top_key = pq.empty() ? std::numeric_limits<float>::min() : pq.topKey();
      } else {
        // gain was updated by success_func in this case
        pq.adjustKey(node, gain_in_pq);
        gpq.top_key = pq.topKey();
      }
      gpq.lock.unlock();
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
        if (lockedModifyPQ(best_id)) return true;
        // if you got a PQ but it fails because the node's gain was wrong or the node couldn't be locked
        // (success_func failed) then we still want to use the standard method
        i = 0;
      }

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
        if (!_pqs[best_id].lock.tryLock()) continue;
        if (lockedModifyPQ(best_id)) return true;
      }
    }

    bool findNextMove() {
      return tryPop();
    }

    Move next_move;

    PartitionedHypergraph& _phg;
    GainCache& _gain_cache;
    const Context& _context;

    vec<rebalancer::GuardedPQ>& _pqs;
    vec<PartitionID>& _target_part;
    vec<rebalancer::NodeState>& _node_state;
    AccessToken _token;
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




/*
 *  TODO parameter tuning experiments
 *  1) test gain vs gain/weight
 *  2) test num PQs to probe
 *  3) scale up num PQs to probe if
 *      a) only little work is left
 *      b) many threads are active
 *  4) test impact of initial distribution scheme and power of k choices for k > 2
 *  5) slow down once we're "close" to being balanced
 */

  template <typename TypeTraits, typename GainTypes>
  void RebalancerV2<TypeTraits, GainTypes>::insertNodesInOverloadedBlocks(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
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
    static constexpr size_t NUM_GAIN_BUCKETS = 12;
    vec<vec<size_t>> freq(num_pqs, vec<size_t>(NUM_GAIN_BUCKETS, 0));

    phg.doParallelForAllNodes([&](HypernodeID u) {
      const PartitionID b = phg.partID(u);
      if (!_is_overloaded[b]) return;

      _node_state[u].markAsMovable();
      auto [target, gain] = impl::computeBestTargetBlock(phg, _context, _gain_cache, u, phg.partID(u));
      _target_part[u] = target;

      size_t bucket_id;
      if (gain > 0.0) bucket_id = 0;
      else if (gain == 0.0) bucket_id = 1;
      else {
        float x = std::log2(-gain);
        bucket_id = std::min<size_t>(2 + x, NUM_GAIN_BUCKETS - 1);
      }


      auto& token = ets_tokens.local();
      int my_pq_id = -1;
      while (true) {
        auto two_ids = token.getTwoRandomPQs();
        my_pq_id = two_ids[0];
        if (freq[two_ids[0]][bucket_id] > freq[two_ids[1]][bucket_id]) {
          my_pq_id = two_ids[1];
        }
        if (_pqs[my_pq_id].lock.tryLock()) {
          break;
        }
      }

      _pqs[my_pq_id].pq.insert(u, gain);
      _pqs[my_pq_id].lock.unlock();
      __atomic_fetch_add(&freq[my_pq_id][bucket_id], 1, __ATOMIC_RELAXED);
      _pq_id[u] = my_pq_id;
    });


    for (rebalancer::GuardedPQ& gpq : _pqs) {
      if (!gpq.pq.empty()) {
        gpq.top_key = gpq.pq.topKey();
      }
    }
  }

  template <typename TypeTraits, typename GainTypes>
  std::pair<int64_t, size_t> RebalancerV2<TypeTraits, GainTypes>::findMoves(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);
    int64_t attributed_gain = 0;
    size_t global_move_id = 0;
    size_t num_overloaded_blocks = _overloaded_blocks.size();

    auto task = [&](size_t task_id) {
      vec<HyperedgeID> edges_with_gain_changes;
      Gain local_attributed_gain = 0;
      vec<vec<HypernodeID>> nodes_to_update(_pqs.size());
      vec<int> pqs_to_update;

      const int seed = phg.initialNumNodes() + task_id;

      impl::NextMoveFinder next_move_finder(seed, _context, phg, _gain_cache, _pqs, _target_part, _node_state);

      while (num_overloaded_blocks > 0 && next_move_finder.findNextMove()) {
        const Move& m = next_move_finder.next_move;
        const PartitionID from = phg.partID(m.node);
        _node_state[m.node].markAsMovedAndUnlock();

        if (phg.partWeight(from) <= _max_part_weights[from]) {
          impl::deactivateOverloadedBlock(&_is_overloaded[from], &num_overloaded_blocks);
          continue;
        }

        edges_with_gain_changes.clear();
        size_t move_id = 0;
        bool moved = phg.changeNodePart(
                      _gain_cache, m.node, m.from, m.to,
                      _max_part_weights[m.to],
                      [&] { move_id = __atomic_fetch_add(&global_move_id, 1, __ATOMIC_RELAXED); },
                      [&](const SyncronizedEdgeUpdate& sync_update) {
                        local_attributed_gain += AttributedGains::gain(sync_update);
                        if (!phg.is_graph && GainCache::triggersDeltaGainUpdate(sync_update)) {
                          edges_with_gain_changes.push_back(sync_update.he);
                        }
                      }
                    );



        if (!moved) continue;


        auto update_neighbor = [&](HypernodeID v) {
          if (v != m.node && _node_state[v].tryLock()) {
            int my_pq_id = _pq_id[v];
            assert(my_pq_id != -1);
            if (nodes_to_update[my_pq_id].empty()) {
              pqs_to_update.push_back(my_pq_id);
            }
            nodes_to_update[my_pq_id].push_back(v);
            next_move_finder.recomputeTopGainMove(v, m);
          }
        };

        // update neighbors
        if constexpr (phg.is_graph) {
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
                  Gain new_gain_int = _gain_cache.gain(v, phg.partID(v), _target_part[v]);
                  float new_gain = impl::transformGain(new_gain_int, phg.nodeWeight(v));
                  pq.adjustKey(v, new_gain);
                }
                // TODO when should we unlock nodes? right after its target part is set? or only now?
                // only now means it can't be moved when it's about to receive an update and
                // it can't be updated when it's about to be moved.
                _node_state[v].unlock();
              }

              gpq.lock.unlock();
              pqs_to_update[i] = pqs_to_update.back();
              pqs_to_update.pop_back();
              nodes_to_update[my_pq_id].clear();
            }
          }
        }

        _moves[move_id] = m;
      }
      __atomic_fetch_add(&attributed_gain, local_attributed_gain, __ATOMIC_RELAXED);
    };

    tbb::task_group tg;
    for (size_t i = 0; i < _context.shared_memory.num_threads; ++i) { tg.run(std::bind(task, i)); }
    tg.wait();

    return std::make_pair(attributed_gain, global_move_id);
  }

  template <typename TypeTraits, typename GainTypes>
  bool RebalancerV2<TypeTraits, GainTypes>::refineInternalParallel(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                                   vec<vec<Move>>* moves_by_part,
                                                                   Metrics& best_metric) {
    auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

    if (_max_part_weights == nullptr) {
      _max_part_weights = &_context.partition.max_part_weights[0];
    }
    if (!_gain_cache.isInitialized()) {
      _gain_cache.initializeGainCache(phg);
    }

    _overloaded_blocks.clear();
    _is_overloaded.assign(phg.k(), false);
    for (PartitionID k = 0; k < phg.k(); ++k) {
      if (phg.partWeight(k) > _max_part_weights[k]) {
        _overloaded_blocks.push_back(k);
        _is_overloaded[k] = 1;
      }
    }

    insertNodesInOverloadedBlocks(hypergraph);

    auto [attributed_gain, num_moves_performed] = findMoves(hypergraph);

    if (moves_by_part != nullptr) {
      moves_by_part->resize(phg.k());
      for (auto& direction : *moves_by_part) direction.clear();
      for (size_t i = 0; i < num_moves_performed; ++i) {
        (*moves_by_part)[_moves[i].from].push_back(_moves[i]);
      }
    }

    best_metric.quality += attributed_gain;
    best_metric.imbalance = metrics::imbalance(phg, _context);

    size_t num_overloaded_blocks = 0;
    for (PartitionID b = 0; b < phg.k(); ++b) {
      if (phg.partWeight(b) > _max_part_weights[b]) {
        num_overloaded_blocks++;
      }
    }

    _max_part_weights = nullptr;

    phg.doParallelForAllNodes([&](HypernodeID u) {
      _node_state[u].reset();
    });

    for (auto& gpq : _pqs) {
      gpq.pq.clear();
    }

    return num_overloaded_blocks == 0;
  }


template <typename TypeTraits, typename GainTypes>
RebalancerV2<TypeTraits, GainTypes>::RebalancerV2(
        HypernodeID num_nodes, const Context& context, GainCache& gain_cache) :
        _context(context),
        _max_part_weights(nullptr),
        _gain_cache(gain_cache),
        _current_k(_context.partition.k),
        _gain(context),
        _moves(num_nodes),
        _target_part(num_nodes, kInvalidPartition),
        _pq_handles(num_nodes, invalid_position),
        _pq_id(num_nodes, -1),
        _node_state(num_nodes) { }

template <typename TypeTraits, typename GainTypes>
RebalancerV2<TypeTraits, GainTypes>::RebalancerV2(
        HypernodeID num_nodes, const Context& context, gain_cache_t gain_cache) :
        RebalancerV2(num_nodes, context, GainCachePtr::cast<GainCache>(gain_cache)) { }


template <typename TypeTraits, typename GainTypes>
bool RebalancerV2<TypeTraits, GainTypes>::refineImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                const vec<HypernodeID>& , Metrics& best_metrics, double) {
  return refineInternalParallel(hypergraph, nullptr, best_metrics);
}

template <typename TypeTraits, typename GainTypes>
void RebalancerV2<TypeTraits, GainTypes>::initializeImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph) {
  auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  unused(phg);
}


template <typename TypeTraits, typename GainTypes>
void RebalancerV2<TypeTraits, GainTypes>::setMaxPartWeightsForRoundImpl(const std::vector<HypernodeWeight>& max_part_weights) {
  _max_part_weights = &max_part_weights[0];
}

template <typename TypeTraits, typename GainTypes>
bool RebalancerV2<TypeTraits, GainTypes>::refineAndOutputMovesImpl(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                              const vec<HypernodeID>& ,
                              vec<vec<Move>>& moves_by_part,
                              Metrics& best_metrics,
                              const double) {
  return refineInternalParallel(hypergraph, &moves_by_part, best_metrics);
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
  #define REBALANCER_V2(X, Y) RebalancerV2<X, Y>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(REBALANCER_V2)

}   // namespace mt_kahypar

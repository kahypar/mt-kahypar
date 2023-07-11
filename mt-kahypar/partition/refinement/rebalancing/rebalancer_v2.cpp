#include <mt-kahypar/datastructures/priority_queue.h>
#include <optional>
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_v2.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/partition/context.h"

#include "pcg_random.hpp"

namespace mt_kahypar {

namespace {
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
    for (PartitionID i = 0; i < context.partition.k; ++i) {
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
} // anonymous namespace


class PQLayout {
public:
  PQLayout(const Context& context, size_t num_nodes) :
          context(context),
          target_part(num_nodes, kInvalidPartition),
          pq_handles(num_nodes, invalid_position),
          pqs(static_cast<size_t>(context.partition.k), VertexPriorityQueue(pq_handles.data(), num_nodes)) { }

  template<typename PartitionedHypergraph, typename GainCache>
  void insertIntoPQ(const PartitionedHypergraph& phg,
                    const GainCache& gain_cache,
                    const HypernodeID v) {
    const PartitionID pv = phg.partID(v);
    ASSERT(pv < context.partition.k);
    auto [target, gain] = computeBestTargetBlock(phg, context, gain_cache, v, pv);
    ASSERT(target < context.partition.k);
    target_part[v] = target;
    pqs[pv].insert(v, gain);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  void updateGain(const PartitionedHypergraph& phg,
                  const GainCache& gain_cache,
                  const HypernodeID v,
                  const Move& move) {
    const PartitionID pv = phg.partID(v);
    ASSERT(pqs[pv].contains(v));
    const PartitionID designatedTargetV = target_part[v];
    float gain = 0;
    PartitionID newTarget = kInvalidPartition;

    if (context.partition.k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, context, gain_cache, v, pv);
    } else {
      std::tie(newTarget, gain) = bestOfThree(phg, context, gain_cache,
                                              v, pv, { designatedTargetV, move.from, move.to });
    }
    target_part[v] = newTarget;
    pqs[pv].adjustKey(v, gain);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  bool findNextMove(const PartitionedHypergraph& phg,
                    const GainCache& gain_cache,
                    Move& m, PartitionID from) {

    if (pqs[from].empty()) return false;

    while (true) {
      const HypernodeID u = pqs[from].top();
      const float estimated_gain = pqs[from].topKey();
      auto [to, gain] = computeBestTargetBlock(phg, context, gain_cache, u, phg.partID(u));

      if (gain >= estimated_gain) {
        m.node = u; m.to = to; m.from = from;
        m.gain = gain;
        pqs[from].deleteTop();
        return true;
      } else {
        target_part[u] = to;
        pqs[from].adjustKey(u, gain);
      }
    }
  }

  PartitionID findNextBlockToMoveFrom() {
    PartitionID best_block = -1;
    float max_gain = std::numeric_limits<float>::min();
    for (PartitionID k : overloaded_blocks) {
      if (!pqs[k].empty() && pqs[k].topKey() > max_gain) {
        max_gain = pqs[k].topKey();
        best_block = k;
      }
    }
    return best_block;
  }

  const Context& context;
  vec<PartitionID> overloaded_blocks;
  vec<PartitionID> target_part;
  vec<PosT> pq_handles;
  using VertexPriorityQueue = ds::MaxHeap<float, HypernodeID>;
  vec<VertexPriorityQueue> pqs;
};



template <typename TypeTraits, typename GainTypes>
RebalancerV2<TypeTraits, GainTypes>::RebalancerV2(const Context& context,
                      GainCache& gain_cache) :
        _context(context),
        _max_part_weights(nullptr),
        _gain_cache(gain_cache),
        _current_k(_context.partition.k),
        _gain(context)
{ }

template <typename TypeTraits, typename GainTypes>
RebalancerV2<TypeTraits, GainTypes>::RebalancerV2(const Context& context,
                      gain_cache_t gain_cache) :
        RebalancerV2(context, GainCachePtr::cast<GainCache>(gain_cache)) { }


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



template <typename TypeTraits, typename GainTypes>
bool RebalancerV2<TypeTraits, GainTypes>::refineInternal(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                    vec<vec<Move>>* moves_by_part,
                    Metrics& best_metric) {
  auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);
  unused(phg);
  unused(moves_by_part); unused(best_metric);
  vec<Move> moves;
  Gain attributed_gain = 0;
  vec<HyperedgeID> edges_with_gain_changes;

  PQLayout pq_layout(_context, phg.initialNumNodes());
  vec<HypernodeID> neighborDeduplicator(phg.initialNumNodes(), 0);
  HypernodeID deduplicationTime = 1;

  for (PartitionID k = 0; k < phg.k(); ++k) {
    if (phg.partWeight(k) > _max_part_weights[k]) {
      pq_layout.overloaded_blocks.push_back(k);
    }
  }

  bool select_block_by_gain = true;

  if (!select_block_by_gain) {
    // work on lighter blocks first (placed at the end)
    std::sort(pq_layout.overloaded_blocks.begin(), pq_layout.overloaded_blocks.end(), [&](PartitionID l, PartitionID r) {
      return phg.partWeight(l) > phg.partWeight(r);
    });
  }

  {
    std::vector<bool> is_overloaded(phg.k(), false);
    for (PartitionID k : pq_layout.overloaded_blocks) is_overloaded[k] = true;
    for (HypernodeID u : phg.nodes()) {
      if (is_overloaded[phg.partID(u)]) {
        pq_layout.insertIntoPQ(phg, _gain_cache, u);
      }
    }
  }

  while (!pq_layout.overloaded_blocks.empty()) {
    PartitionID from;
    if (select_block_by_gain) {
      from = pq_layout.findNextBlockToMoveFrom();
      if (from == -1) break;
    } else {
      from = pq_layout.overloaded_blocks.back();
    }


    Move m;
    if (!pq_layout.findNextMove(phg, _gain_cache, m, from)) {
      break;
      // TODO this shouldn't happen
    }

    edges_with_gain_changes.clear();
    phg.changeNodePart(
            _gain_cache, m.node, m.from, m.to,
            std::numeric_limits<HypernodeWeight>::max(),
            [] { },
            [&](HyperedgeID he, HyperedgeWeight edge_weight, HypernodeID edge_size,
                HypernodeID pin_count_in_from_part_after, HypernodeID pin_count_in_to_part_after) {
              attributed_gain += AttributedGains::gain(he, edge_weight, edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after);
              if (GainCache::triggersDeltaGainUpdate(edge_size, pin_count_in_from_part_after, pin_count_in_to_part_after)) {
                edges_with_gain_changes.push_back(he);
              }
            });

    moves.push_back(m);

    for (HyperedgeID e : edges_with_gain_changes) {
      if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          // TODO beware. this can reinsert nodes
          if (pq_layout.pq_handles[v] != invalid_position && neighborDeduplicator[v] != deduplicationTime) {
            pq_layout.updateGain(phg, _gain_cache, v, m);
            neighborDeduplicator[v] = deduplicationTime;  // not needed for graphs... implement specialized loop?
          }
        }
      }
    }

    if (++deduplicationTime == 0) {
      neighborDeduplicator.assign(neighborDeduplicator.size(), 0);
      deduplicationTime = 1;
    }

    if (phg.partWeight(from) <= _max_part_weights[from]) {
      auto it = std::find(pq_layout.overloaded_blocks.begin(), pq_layout.overloaded_blocks.end(), from);
      pq_layout.overloaded_blocks.erase(it);
    }
    if (phg.partWeight(m.to) > _max_part_weights[m.to]) {
      pq_layout.overloaded_blocks.push_back(m.to);
      if (!select_block_by_gain) {
        std::sort(pq_layout.overloaded_blocks.begin(), pq_layout.overloaded_blocks.end(), [&](PartitionID l, PartitionID r) {
          return phg.partWeight(l) > phg.partWeight(r);
        });
      }
    }
  }

  if (moves_by_part != nullptr) {
    moves_by_part->resize(phg.k());
    for (auto& direction : *moves_by_part) direction.clear();
    for (const Move& m : moves) {
      (*moves_by_part)[m.from].push_back(m);
    }
  }

  best_metric.quality += attributed_gain;
  best_metric.imbalance = metrics::imbalance(phg, _context);

  return pq_layout.overloaded_blocks.empty();
}

namespace impl {

// TODO can be moved to inner member of ParallelPQ
  struct GuardedPQ {
    GuardedPQ(PosT *handles, size_t num_nodes) : pq(handles, num_nodes) { }

    SpinLock lock;
    ds::MaxHeap<float, HypernodeID> pq;
    float top_key = std::numeric_limits<float>::min();
  };

  struct AccessToken {
    AccessToken(int seed, size_t num_pqs) : dist(0, num_pqs - 1) {
      rng.seed(seed);
    }

    size_t getRandomPQ() {
      return dist(rng);
    }

    std::array<size_t, 2> getTwoRandomPQs() {
      std::array<size_t, 2> result({getRandomPQ(), getRandomPQ()});
      while (result[0] != result[1]) { result[1] = getRandomPQ(); }
      return result;
    }

    pcg32 rng;
    std::uniform_int_distribution<size_t> dist;
  };

  struct NodeState {
    uint8_t state = 0;

    bool canMove() const { return state == 1; }

    bool isLocked() const { return state == 2; }

    bool wasMoved() const { return state == 3; }

    // Returns true if the node is marked as movable, is not locked and taking the lock now succeeds
    bool tryLock() {
      uint8_t expected = 1;
      return state == 1 && __atomic_compare_exchange_n(&state, &expected, 2, false, __ATOMIC_ACQUIRE, __ATOMIC_RELAXED);
    }

    void unlock() { __atomic_store_n(&state, 1, __ATOMIC_RELEASE); }

    void markAsMovedAndUnlock() { __atomic_store_n(&state, 3, __ATOMIC_RELEASE); }

    void markAsMovable() { state = 1; }
  };

  struct ParallelPQ {
    // only alloc PQs if the block is overloaded
    void initialize(size_t num_pqs, PosT *handles, size_t num_nodes) {
      pqs.assign(num_pqs, GuardedPQ(handles, num_nodes));
    }

    size_t push(AccessToken& token, HypernodeID node, float gain) {
      size_t pq_id = token.getRandomPQ();
      while (!pqs[pq_id].lock.tryLock()) {
        pq_id = token.getRandomPQ();
      }
      pqs[pq_id].pq.insert(node, gain);
      pqs[pq_id].top_key = pqs[pq_id].pq.topKey();
      pqs[pq_id].lock.unlock();
      return pq_id;
    }

    template<typename SuccessFunc>
    bool lockedModifyPQ(size_t best_id, SuccessFunc success_func) {
      HypernodeID node = pqs[best_id].pq.top();
      float gain = pqs[best_id].pq.topKey();
      bool success = success_func(node, gain);
      if (success) {
        pqs[best_id].pq.deleteTop();
        if (!pqs[best_id].pq.empty()) {
          pqs[best_id].top_key = pqs[best_id].pq.topKey();
        } else {
          pqs[best_id].top_key = std::numeric_limits<float>::min();
        }
      } else {
        // gain was updated by success_func in this case
        pqs[best_id].pq.adjustKey(node, gain);
        pqs[best_id].top_key = pqs[best_id].pq.topKey();
      }
      pqs[best_id].lock.unlock();
      return success;
    }

    template<typename SuccessFunc>
    bool tryPop(AccessToken& token, SuccessFunc success_func) {
      static constexpr size_t num_tries = 32;
      for (size_t i = 0; i < num_tries; ++i) {
        auto two = token.getTwoRandomPQs();
        auto& first = pqs[two[0]];
        auto& second = pqs[two[1]];
        if (first.pq.empty() && second.pq.empty()) continue;
        size_t best_id = two[0];
        if (first.pq.empty() || first.top_key < second.top_key) best_id = two[1];
        if (!pqs[best_id].lock.tryLock()) continue;
        // could also check for top key. would want to distinguish tries that failed due to high contention
        // vs approaching the end
        if (pqs[best_id].pq.empty()) {
          pqs[best_id].lock.unlock();
          continue;
        }
        if (lockedModifyPQ(best_id, success_func)) return true;
        // if you got a PQ but it fails because the node's gain was wrong or the node couldn't be locked
        // (success_func failed) then we still want to use the standard method
        i = 0;
      }

      while (true) {
        float best_key = std::numeric_limits<float>::min();
        int best_id = -1;
        for (size_t i = 0; i < pqs.size(); ++i) {
          if (!pqs[i].pq.empty() && pqs[i].top_key > best_key) {
            best_key = pqs[i].top_key;
            best_id = i;
          }
        }
        if (best_id == -1) return false;
        if (!pqs[best_id].lock.tryLock()) continue;
        if (lockedModifyPQ(best_id, success_func)) return true;
      }
    }

    vec <GuardedPQ> pqs;
  };


// TODO move all this to members of the rebalancer
  class ParallelPQLayout {
  public:
    ParallelPQLayout(const Context& context, size_t num_nodes) :
            context(context),
            target_part(num_nodes, kInvalidPartition),
            pq_handles(num_nodes, invalid_position),
            pq_id(num_nodes, -1),
            node_state(num_nodes) { }

    void initializeRelaxedPQs(size_t num_pqs) {
      pq.initialize(num_pqs, pq_handles.data(), pq_handles.size());
    }

    template<typename PartitionedHypergraph, typename GainCache>
    void recomputeTopGainMove(const PartitionedHypergraph& phg, const GainCache& gain_cache, HypernodeID v,
                              const Move& move) {
      const PartitionID designatedTargetV = target_part[v];
      float gain = 0;
      PartitionID newTarget = kInvalidPartition;
      if (context.partition.k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
        std::tie(newTarget, gain) = computeBestTargetBlock(phg, context, gain_cache, v, phg.partID(v));
      } else {
        std::tie(newTarget, gain) = bestOfThree(phg, context, gain_cache,
                                                v, phg.partID(v), {designatedTargetV, move.from, move.to});
      }
      target_part[v] = newTarget;
    }


    template<typename PartitionedHypergraph, typename GainCache>
    bool findNextMove(const PartitionedHypergraph& phg,
                      const GainCache& gain_cache,
                      Move& m, AccessToken& token) {
      return pq.tryPop(token, [&](HypernodeID u, float& gain_in_pq) -> bool {
        if (!node_state[u].tryLock()) return false;
        auto[to, true_gain] = computeBestTargetBlock(phg, context, gain_cache, u, phg.partID(u));
        if (true_gain >= gain_in_pq) {
          m.node = u;
          m.to = to;
          m.from = phg.partID(u);
          m.gain = true_gain;
          return true;
        } else {
          target_part[u] = to;
          // the tryPop function will use this value to adjust the key.
          gain_in_pq = true_gain;
          node_state[u].unlock();
          return false;
        }
      });
    }

    const Context& context;
    vec <PartitionID> overloaded_blocks;
    vec <uint8_t> is_overloaded;
    vec <PartitionID> target_part;
    vec <PosT> pq_handles;
    vec<int> pq_id;
    vec <NodeState> node_state;
    ParallelPQ pq;
  };

}   // namespace impl

template <typename TypeTraits, typename GainTypes>
bool RebalancerV2<TypeTraits, GainTypes>::refineInternalParallel(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                         vec<vec<Move>>* moves_by_part,
                                                         Metrics& best_metric) {
  auto& phg = utils::cast<PartitionedHypergraph>(hypergraph);

  const auto& max_part_weights = _context.partition.max_part_weights;
  if (_max_part_weights == nullptr) {
    _max_part_weights = &_context.partition.max_part_weights[0];
  }
  if (!_gain_cache.isInitialized()) {
    _gain_cache.initializeGainCache(phg);
  }

  vec<Move> moves(phg.initialNumNodes());
  Gain attributed_gain = 0;

  impl::ParallelPQLayout pq_layout(_context, phg.initialNumNodes());

  pq_layout.is_overloaded.assign(phg.k(), false);
  for (PartitionID k = 0; k < phg.k(); ++k) {
    if (phg.partWeight(k) > _max_part_weights[k]) {
      pq_layout.overloaded_blocks.push_back(k);
      pq_layout.is_overloaded[k] = 1;
    }
  }

  const size_t num_pqs = 2 * _context.shared_memory.num_threads;
  pq_layout.initializeRelaxedPQs(num_pqs);

  std::atomic<int> seed { 555 };
  tbb::enumerable_thread_specific<impl::AccessToken> ets_tokens([&]() {
    return impl::AccessToken(seed.fetch_add(1, std::memory_order_relaxed), num_pqs);
  });
  static constexpr size_t NUM_GAIN_BUCKETS = 12;
  vec<vec<size_t>> freq(num_pqs, vec<size_t>(NUM_GAIN_BUCKETS, 0));

  phg.doParallelForAllNodes([&](HypernodeID u) {
    PartitionID b = phg.partID(u);

    if (!pq_layout.is_overloaded[b]) return;

    pq_layout.node_state[u].markAsMovable();
    auto [target, gain] = computeBestTargetBlock(phg, _context, _gain_cache, u, phg.partID(u));
    pq_layout.target_part[u] = target;

    size_t bucket_id;
    if (gain > 0.0) bucket_id = 0;
    else if (gain == 0.0) bucket_id = 1;
    else {
      float x = std::log2(-gain);
      bucket_id = std::min<size_t>(2 + x, NUM_GAIN_BUCKETS - 1);
    }


    auto& token = ets_tokens.local();
    int pq_id = -1;
    while (true) {
      auto two_ids = token.getTwoRandomPQs();
      pq_id = two_ids[0];
      if (freq[two_ids[0]][bucket_id] > freq[two_ids[1]][bucket_id]) {
        pq_id = two_ids[1];
      }
      if (pq_layout.pq.pqs[pq_id].lock.tryLock()) {
        break;
      }
    }

    pq_layout.pq.pqs[pq_id].pq.insert(u, gain);
    pq_layout.pq.pqs[pq_id].lock.unlock();
    __atomic_fetch_add(&freq[pq_id][bucket_id], 1, __ATOMIC_RELAXED);
    pq_layout.pq_id[u] = pq_id;
  });


  for (size_t pq_id = 0; pq_id < num_pqs; ++pq_id) {
    if (!pq_layout.pq.pqs[pq_id].pq.empty()) {
      pq_layout.pq.pqs[pq_id].top_key = pq_layout.pq.pqs[pq_id].pq.topKey();
    }
  }


  size_t global_move_id = 0;
  size_t num_overloaded_blocks = pq_layout.overloaded_blocks.size();
  auto task = [&](size_t ) {
    vec<HyperedgeID> edges_with_gain_changes;
    Gain local_attributed_gain = 0;
    vec<vec<HypernodeID>> nodes_to_update(num_pqs);
    vec<int> pqs_to_update;
    auto& token = ets_tokens.local();
    Move m;
    while (num_overloaded_blocks > 0 && pq_layout.findNextMove(phg, _gain_cache, m, token)) {
      const PartitionID from = phg.partID(m.node);
      pq_layout.node_state[m.node].markAsMovedAndUnlock();
      if (phg.partWeight(from) <= max_part_weights[from]) {
        if (pq_layout.is_overloaded[from]) {
          uint8_t expected = 1;
          if (__atomic_compare_exchange_n(
                  &pq_layout.is_overloaded[from], &expected, 0, false, __ATOMIC_ACQUIRE, __ATOMIC_RELAXED)) {
            __atomic_fetch_sub(&num_overloaded_blocks, 1, __ATOMIC_RELAXED);
          }
        }
        continue;
      }

      edges_with_gain_changes.clear();
      bool moved = phg.changeNodePart(
              _gain_cache, m.node, m.from, m.to,
              _max_part_weights[m.to],
              [] { },
              [&](HyperedgeID he, HyperedgeWeight edge_weight, HypernodeID edge_size,
                  HypernodeID pin_count_in_from_part_after, HypernodeID pin_count_in_to_part_after) {
                local_attributed_gain += AttributedGains::gain(he, edge_weight, edge_size, pin_count_in_from_part_after,
                                                               pin_count_in_to_part_after);
                if (GainCache::triggersDeltaGainUpdate(edge_size, pin_count_in_from_part_after,
                                                       pin_count_in_to_part_after)) {
                  edges_with_gain_changes.push_back(he);
                }
              }
      );

      if (!moved) continue;

      size_t move_id = __atomic_fetch_add(&global_move_id, 1, __ATOMIC_RELAXED);
      moves[move_id] = m;

      // TODO write specialized loop for graphs
      for (HyperedgeID e : edges_with_gain_changes) {
        if (phg.edgeSize(e) < _context.partition.ignore_hyperedge_size_threshold) {
          for (HypernodeID v : phg.pins(e)) {
            if (v != m.node && pq_layout.node_state[v].tryLock()) {
              int pq_id = pq_layout.pq_id[v];
              assert(pq_id != -1);
              if (nodes_to_update[pq_id].empty()) {
                pqs_to_update.push_back(pq_id);
              }
              nodes_to_update[pq_id].push_back(v);
              pq_layout.recomputeTopGainMove(phg, _gain_cache, v, m);
            }
          }
        }
      }


      while (!pqs_to_update.empty()) {
        for (size_t i = 0; i < pqs_to_update.size(); ++i) {
          int pq_id = pqs_to_update[i];
          auto& gpq = pq_layout.pq.pqs[pq_id];
          auto& pq = gpq.pq;
          if (gpq.lock.tryLock()) {
            for (HypernodeID v : nodes_to_update[pq_id]) {
              if (pq.contains(v)) {
                Gain new_gain_int = _gain_cache.gain(v, phg.partID(v), pq_layout.target_part[v]);
                float new_gain = transformGain(new_gain_int, phg.nodeWeight(v));
                pq.adjustKey(v, new_gain);
              }
              // TODO when should we unlock nodes? right after its target part is set? or only now?
              // only now means it can't be moved when it's about to receive an update and
              // it can't be updated when it's about to be moved.
              pq_layout.node_state[v].unlock();
            }

            gpq.lock.unlock();
            pqs_to_update[i] = pqs_to_update.back();
            pqs_to_update.pop_back();
            nodes_to_update[pq_id].clear();
          }
        }
      }

      __atomic_fetch_add(&attributed_gain, local_attributed_gain, __ATOMIC_RELAXED);
    }
  };

  tbb::task_group tg;
  for (size_t i = 0; i < _context.shared_memory.num_threads; ++i) { tg.run(std::bind(task, i)); }
  tg.wait();

  if (moves_by_part != nullptr) {
    moves.resize(global_move_id);
    moves_by_part->resize(phg.k());
    for (auto& direction : *moves_by_part) direction.clear();
    for (const Move& m : moves) {
      (*moves_by_part)[m.from].push_back(m);
    }
  }

  best_metric.quality += attributed_gain;
  best_metric.imbalance = metrics::imbalance(phg, _context);

  _max_part_weights = nullptr;

  return num_overloaded_blocks == 0;
}


// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
  #define REBALANCER_V2(X, Y) RebalancerV2<X, Y>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(REBALANCER_V2)

}   // namespace mt_kahypar

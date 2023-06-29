#include <mt-kahypar/datastructures/priority_queue.h>
#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_v2.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

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
    auto [target, gain] = computeBestTargetBlock(phg, gain_cache, v, pv);
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
    Gain gain = 0;
    PartitionID newTarget = kInvalidPartition;

    if (context.partition.k < 4 || designatedTargetV == move.from || designatedTargetV == move.to) {
      std::tie(newTarget, gain) = computeBestTargetBlock(phg, gain_cache, v, pv);
    } else {
      std::tie(newTarget, gain) = bestOfThree(phg, gain_cache,
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
      const Gain estimated_gain = pqs[from].topKey();
      auto [to, gain] = computeBestTargetBlock(phg, gain_cache, u, phg.partID(u));

      if (gain >= estimated_gain) {
        m.node = u; m.to = to; m.from = from;
        m.gain = gain;
        pqs[from].deleteTop();
        return true;
      } else {
        pqs[from].adjustKey(u, gain);
      }
    }
  }

  PartitionID findNextBlockToMoveFrom() {
    PartitionID best_block = -1;
    Gain max_gain = std::numeric_limits<Gain>::min();
    for (PartitionID k : overloaded_blocks) {
      if (!pqs[k].empty() && pqs[k].topKey() > max_gain) {
        max_gain = pqs[k].topKey();
        best_block = k;
      }
    }
    return best_block;
  }

  template<typename PartitionedHypergraph, typename GainCache>
  void deltaGainUpdates(PartitionedHypergraph& phg,
                        GainCache& gain_cache,
                        const HyperedgeID he,
                        const HyperedgeWeight edge_weight,
                        const PartitionID from,
                        const HypernodeID pin_count_in_from_part_after,
                        const PartitionID to,
                        const HypernodeID pin_count_in_to_part_after) {
    gain_cache.deltaGainUpdate(phg, he, edge_weight, from,
                               pin_count_in_from_part_after, to, pin_count_in_to_part_after);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  std::pair<PartitionID, HyperedgeWeight> computeBestTargetBlock(const PartitionedHypergraph& phg,
                                                                 const GainCache& gain_cache,
                                                                 const HypernodeID u,
                                                                 const PartitionID from) {
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
      if (gain > 0 && multiply_positive_gains) {
        gain *= wu;
      } else if (gain < 0) {
        gain /= wu;
      }
    }

    return std::make_pair(to, gain);
  }

  template<typename PartitionedHypergraph, typename GainCache>
  std::pair<PartitionID, HyperedgeWeight> bestOfThree(const PartitionedHypergraph& phg,
                                                      const GainCache& gain_cache,
                                                      HypernodeID u,
                                                      PartitionID from,
                                                      std::array<PartitionID, 3> parts) {

    const HypernodeWeight wu = phg.nodeWeight(u);
    const HypernodeWeight from_weight = phg.partWeight(from);
    PartitionID to = kInvalidPartition;
    HyperedgeWeight to_benefit = std::numeric_limits<HyperedgeWeight>::min();
    HypernodeWeight best_to_weight = from_weight - wu;
    for (PartitionID i : parts) {
      if (i != from && i != kInvalidPartition) {
        const HypernodeWeight to_weight = phg.partWeight(i);
        const HyperedgeWeight penalty = gain_cache.benefitTerm(u, i);
        if ( ( penalty > to_benefit || ( penalty == to_benefit && to_weight < best_to_weight ) ) &&
             to_weight + wu <= context.partition.max_part_weights[i] ) {
          to_benefit = penalty;
          to = i;
          best_to_weight = to_weight;
        }
      }
    }

    Gain gain = std::numeric_limits<Gain>::min();
    if (to != kInvalidPartition) {
      gain = to_benefit - gain_cache.penaltyTerm(u, phg.partID(u));
      if (gain > 0 && multiply_positive_gains) {
        gain *= wu;
      } else if (gain < 0) {
        gain /= wu;
      }
    }

    return std::make_pair(to, gain);
  }

  const Context& context;
  vec<PartitionID> overloaded_blocks;
  vec<PartitionID> target_part;
  vec<PosT> pq_handles;
  using VertexPriorityQueue = ds::MaxHeap<Gain, HypernodeID>;
  vec<VertexPriorityQueue> pqs;
  bool multiply_positive_gains = true;
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
                const vec<HypernodeID>& ,
                Metrics& best_metrics,
                double) {
  return refineInternal(hypergraph, nullptr, best_metrics);
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
  return refineInternal(hypergraph, &moves_by_part, best_metrics);
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
      // phg.changeNodePart(_gain_cache, m.node, m.to, m.from);
      // TODO interestingly we assume that all of these moves are applied to the partition still
      // --> don't revert. the interleaveMoveSequence function even applies the last unneeded rebalancing moves
      // and lets the rollback take them back again
    }
  }

  best_metric.quality += attributed_gain;
  best_metric.imbalance = metrics::imbalance(phg, _context);

  return pq_layout.overloaded_blocks.empty();
}

template <typename TypeTraits, typename GainTypes>
bool RebalancerV2<TypeTraits, GainTypes>::refineInternalParallel(mt_kahypar_partitioned_hypergraph_t& hypergraph,
                                                         vec<vec<Move>>* moves_by_part,
                                                         Metrics& best_metric) {
  return true;
}


// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
  #define REBALANCER_V2(X, Y) RebalancerV2<X, Y>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(REBALANCER_V2)

}   // namespace mt_kahypar

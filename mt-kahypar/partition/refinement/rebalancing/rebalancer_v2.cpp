#include "mt-kahypar/partition/refinement/rebalancing/rebalancer_v2.h"

#include "mt-kahypar/partition/refinement/gains/gain_definitions.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

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
  return false;
}


// explicitly instantiate so the compiler can generate them when compiling this cpp file
namespace {
  #define REBALANCER_V2(X, Y) RebalancerV2<X, Y>
}

// explicitly instantiate so the compiler can generate them when compiling this cpp file
INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES(REBALANCER_V2)

}   // namespace mt_kahypar

#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/evolutionary/population.h"

namespace mt_kahypar {

// Forward Declaration
class TargetGraph;

template<typename TypeTraits>
class EvoPartitioner {

    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

    public:
    static PartitionedHypergraph partition(Hypergraph& hypergraph,
                                           Context& context,
                                           TargetGraph* target_graph = nullptr);
    static std::string generateInitialPopulation(EvoPartitioner<TypeTraits>::Hypergraph& hg,
                                    Context& context,
                                    TargetGraph* target_graph,
                                    Population& population);
    static const Individual & generateIndividual(EvoPartitioner<TypeTraits>::Hypergraph& hg,
                                            Context& context,
                                            TargetGraph* target_graph,
                                            Population& population);
    static std::string performEvolution(EvoPartitioner<TypeTraits>::Hypergraph& hg,
                                    Context& context,
                                    TargetGraph* target_graph,
                                    Population& population);

    private:
    static EvoDecision decideNextMove(const Context& context);
    static EvoMutateStrategy decideNextMutation(const Context& context);
    static vec<PartitionID> combinePartitions(const Context& context, Population& population, std::vector<size_t> ids);
    static std::string performCombine(EvoPartitioner<TypeTraits>::Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population);
    static std::string performMutation(EvoPartitioner<TypeTraits>::Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population);
};

}  // namespace mt_kahypar
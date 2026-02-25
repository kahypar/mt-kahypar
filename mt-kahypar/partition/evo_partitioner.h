#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/evolutionary/population.h"
#include <mutex>

namespace mt_kahypar {

// Forward Declaration
class TargetGraph;

template<typename TypeTraits>
class EvoPartitioner : public Partitioner<TypeTraits> {

    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    using Partitioner<TypeTraits>::configurePreprocessing;
    using Partitioner<TypeTraits>::setupContext;
    using Partitioner<TypeTraits>::preprocess;

    public:
    static PartitionedHypergraph partition(Hypergraph& hypergraph,
                                           Context& context,
                                           TargetGraph* target_graph = nullptr);
    static std::string generateInitialPopulation(const Hypergraph& hg,
                                    Context& context,
                                    TargetGraph* target_graph,
                                    Population& population);
    static const Individual & generateIndividual(const Hypergraph& hg,
                                            Context& context,
                                            TargetGraph* target_graph,
                                            Population& population);
    static std::string performEvolution(const Hypergraph& hg,
                                    Context& context,
                                    TargetGraph* target_graph,
                                    Population& population);

    private:
    static std::mutex best_tracking_mutex_;
    static HyperedgeWeight global_best_fitness_;
    static std::chrono::milliseconds global_best_time_;
    static std::string diff_matrix_history;
    static std::mutex diff_matrix_history_mutex;

    static EvoDecision decideNextMove(const Context& context);
    static EvoMutateStrategy decideNextMutation(const Context& context);
    static vec<PartitionID> combinePartitions(const Context& context, Population& population, const std::vector<size_t>& ids);
    static std::string performCombine(const Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population);
    static std::string performMutation(const Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population);
    static std::string checkAndLogNewBest(HyperedgeWeight fitness, const std::string& operation_type, std::chrono::milliseconds current_time);
};

}  // namespace mt_kahypar
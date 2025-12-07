#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/evolutionary/population.h"
#include <mutex>
#include <memory>
#include <vector>

namespace mt_kahypar {

// Forward Declaration
class TargetGraph;

struct ContextModifierParameters {
    bool use_random_partitions = false;
    bool use_degree_sorted_partitions = false;
    PartitionID k = 2;
    double epsilon = 0.03;
    bool recursive_bipartitioning = false;
};

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
                                            Population& population,
                                            bool insert_into_population = true);
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
    static std::string iteration_log_history;
    static std::mutex iteration_log_mutex;

    static thread_local std::vector<std::unique_ptr<Individual>> thread_local_temporaries_;
    static const Individual& addThreadLocalTemporary(Individual&& individual);
    static void clearThreadLocalTemporaries();

    static EvoDecision decideNextMove(const Context& context, std::mt19937* rng = nullptr);
    static ContextModifierParameters decideContextModificationParameters(const Context& context);
    static std::vector<PartitionID> createDegreeSortedPartition(const Hypergraph& hypergraph, const Context& context);
    static std::vector<PartitionID> createRandomPartition(const Hypergraph& hypergraph, const Context& context);
    static EvoMutateStrategy decideNextMutation(const Context& context, std::mt19937* rng = nullptr);
    static vec<PartitionID> combinePartitions(const Context& context, Population& population, const std::vector<size_t>& ids);
    static vec<PartitionID> combineModifiedPartitions(const Context& context, std::vector<std::vector<PartitionID>> parent_partitions);
    static Individual performCombine(const Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population);
    static Individual performModifiedCombine(const Hypergraph& hg, const Context&  context, ContextModifierParameters params, TargetGraph* target_graph, Population& population, std::mt19937* rng = nullptr);
    static Context modifyContext(const Context& context, ContextModifierParameters params);
    static Individual performMutation(const Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population, std::mt19937* rng = nullptr);
    static std::string checkAndLogNewBest(HyperedgeWeight fitness, const std::string& operation_type, std::chrono::milliseconds current_time);
    static std::string EvoPartitioner<TypeTraits>::insert_individual_into_population(Individual&& individual, const Context& context, Population& population);
};

}  // namespace mt_kahypar
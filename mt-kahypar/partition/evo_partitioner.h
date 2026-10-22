#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/evolutionary/population.h"
#include "mt-kahypar/partition/evolutionary/stopping_criteria.h"
#include "mt-kahypar/partition/evolutionary/evo_statistics.h"
#include "mt-kahypar/partition/evolutionary/evo_run_log.h"
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

struct CombineResult {
  std::shared_ptr<Individual> individual;
  EvoCombineStrategy strategy;
  HyperedgeWeight parent_best_fitness = 0;
};

struct MutationResult {
  std::shared_ptr<Individual> individual;
  EvoMutateStrategy strategy;
  HyperedgeWeight parent_best_fitness = 0;
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
    static void generateInitialPopulation(const Hypergraph& hg,
                                    Context& context,
                                    TargetGraph* target_graph,
                                    Population& population,
                                    evolutionary::EvoStatistics& stats,
                                    evolutionary::EvoRunLog& log);
    static const std::shared_ptr<Individual> createInsertIndividual(const Hypergraph& hg,
                                            Context& context,
                                            TargetGraph* target_graph,
                                            Population& population);
    static void performEvolution(const Hypergraph& hg,
                                    Context& context,
                                    TargetGraph* target_graph,
                                    Population& population,
                                    evolutionary::EvoStatistics& stats,
                                    evolutionary::EvoRunLog& log);

    static std::vector<PartitionID> createPartition(const Hypergraph& input_hg, Context& context, TargetGraph* target_graph);

    private:
    friend class EvoOperatorTest;
    static thread_local std::vector<std::unique_ptr<Individual>> thread_local_temporaries_;

    static std::shared_ptr<Individual> createIndividual(const Hypergraph& input_hg, Context& context, TargetGraph* target_graph);
    static CombineResult performCombine(const Hypergraph& hg, const Context& context, TargetGraph* target_graph, Population& population, ContextModifierParameters modified_combine_params,std::mt19937* rng = nullptr);
    static MutationResult performMutation(const Hypergraph &hg, const Context &context,
                                          TargetGraph *target_graph, Population &population,
                                          std::mt19937 *rng = nullptr);
    static bool insert_individual_into_population(
        std::shared_ptr<Individual> individual,
        const Context &context,
        Population &population,
        evolutionary::EvoStatistics& stats,
        evolutionary::EvoRunLog& log,
        int iteration,
        evolutionary::EvoOperator op,
        HyperedgeWeight parent_best_fitness);
    static void reportStatistics(const Context& context,
                                 const evolutionary::EvoStatistics& stats,
                                 const Population& population,
                                 int total_iterations);
};
}  // namespace mt_kahypar
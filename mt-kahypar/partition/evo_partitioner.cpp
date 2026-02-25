#include "evo_partitioner.h"
#include "partitioner.cpp"


namespace mt_kahypar {

    template<typename TypeTraits>
    typename EvoPartitioner<TypeTraits>::PartitionedHypergraph EvoPartitioner<TypeTraits>::partition(
        Hypergraph& hypergraph, Context& context, TargetGraph* target_graph) {
        configurePreprocessing(hypergraph, context);
        setupContext(hypergraph, context, target_graph);

        io::printContext(context);
        io::printMemoryPoolConsumption(context);
        io::printInputInformation(context, hypergraph);

    #ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
        bool map_partition_to_target_graph_at_the_end = false;
        if ( context.partition.objective == Objective::steiner_tree &&
            context.mapping.use_two_phase_approach ) {
        map_partition_to_target_graph_at_the_end = true;
        context.partition.objective = Objective::km1;
        context.setupGainPolicy();
        }
    #endif

        // ################## PREPROCESSING ##################
        utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
        timer.start_timer("preprocessing", "Preprocessing");
        DegreeZeroHypernodeRemover<TypeTraits> degree_zero_hn_remover(context);
        LargeHyperedgeRemover<TypeTraits> large_he_remover(context);
        preprocess(hypergraph, context, target_graph);
        sanitize(hypergraph, context, degree_zero_hn_remover, large_he_remover);
        timer.stop_timer("preprocessing");

        Population population;
        generateInitialPopulation(hypergraph, context, target_graph, population);

        performEvolution(hypergraph, context, target_graph, population);

        PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);
        
        std::vector<PartitionID> best = population.individualAt(population.best()).partition();
        partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
            partitioned_hypergraph.setOnlyNodePart(hn, best[hn]);
        });
        LOG << "Before initialization";
        partitioned_hypergraph.initializePartition();
        LOG << "After initialization";
        /*for ( const HypernodeID& hn : partitioned_hypergraph.nodes() ) {
            partitioned_hypergraph.setOnlyNodePart(hn, best[hn]);
        }
        partitioned_hypergraph.initializePartition();*/

        ASSERT([&] {
            bool success = true;
            if ( partitioned_hypergraph.hasFixedVertices() ) {
                for ( const HypernodeID& hn : partitioned_hypergraph.nodes() ) {
                    if ( partitioned_hypergraph.isFixed(hn) &&
                        partitioned_hypergraph.fixedVertexBlock(hn) != partitioned_hypergraph.partID(hn) ) {
                        LOG << "Node" << hn << "is fixed to block" << partitioned_hypergraph.fixedVertexBlock(hn)
                            << ", but is assigned to block" << partitioned_hypergraph.partID(hn);
                        success = false;
                    }
                }
            }
            return success;
        }(), "Some fixed vertices are not assigned to their corresponding block");

        // ################## POSTPROCESSING ##################
        timer.start_timer("postprocessing", "Postprocessing");
        large_he_remover.restoreLargeHyperedges(partitioned_hypergraph);
        degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hypergraph);
        forceFixedVertexAssignment(partitioned_hypergraph, context);
        timer.stop_timer("postprocessing");

    #ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC
        if ( map_partition_to_target_graph_at_the_end ) {
        ASSERT(target_graph);
        context.partition.objective = Objective::steiner_tree;
        timer.start_timer("one_to_one_mapping", "One-To-One Mapping");
        InitialMapping<TypeTraits>::mapToTargetGraph(
            partitioned_hypergraph, *target_graph, context);
        timer.stop_timer("one_to_one_mapping");
        }
    #endif

        if (context.partition.verbose_output) {
            io::printHypergraphInfo(partitioned_hypergraph.hypergraph(), context,
                                    "Uncoarsened Hypergraph", context.partition.show_memory_consumption);
            io::printStripe();
        }

        return partitioned_hypergraph;
    }


    template<typename TypeTraits>
    void EvoPartitioner<TypeTraits>::generateInitialPopulation(EvoPartitioner<TypeTraits>::Hypergraph& hg, Context& context, TargetGraph* target_graph, Population& population) { 
        int timelimit = context.partition.time_limit;
        //context.evolutionary.dynamic_population_size = true;
        //context.evolutionary.population_size = 50;
        utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
        // INITIAL POPULATION
        if (context.evolutionary.dynamic_population_size) {
            HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
            timer.start_timer("evolutionary", "Evolutionary");
            generateIndividual(hg, context, target_graph, population);
            timer.stop_timer("evolutionary");

            ++context.evolutionary.iteration;
            int dynamic_population_size = std::round(context.evolutionary.dynamic_population_amount_of_time
                                                    * context.partition.time_limit
                                                    / timer.get("evolutionary"));
            int minimal_size = std::max(dynamic_population_size, 3);

            context.evolutionary.population_size = std::min(minimal_size, 50);
            LOG << context.evolutionary.population_size;
            LOG << population;
        }
        context.evolutionary.edge_frequency_amount = sqrt(context.evolutionary.population_size);
        LOG << "EDGE-FREQUENCY-AMOUNT";
        LOG << context.evolutionary.edge_frequency_amount;
        while (population.size() < context.evolutionary.population_size &&
            timer.get("evolutionary") <= timelimit) {
            ++context.evolutionary.iteration;
            timer.start_timer("evolutionary", "Evolutionary");
            generateIndividual(hg, context, target_graph, population);
            timer.stop_timer("evolutionary");
            //LOG << population;
        }
    }

    template<typename TypeTraits>
    const Individual & EvoPartitioner<TypeTraits>::generateIndividual(EvoPartitioner<TypeTraits>::Hypergraph& hypergraph, Context& context, TargetGraph* target_graph, Population& population) {
        EvoPartitioner<TypeTraits>::PartitionedHypergraph partitioned_hypergraph;
        if (context.partition.mode == Mode::direct) {
            partitioned_hypergraph = Multilevel<TypeTraits>::partition(hypergraph, context, target_graph);
        } else if (context.partition.mode == Mode::recursive_bipartitioning) {
            partitioned_hypergraph = RecursiveBipartitioning<TypeTraits>::partition(hypergraph, context, target_graph);
        } else if (context.partition.mode == Mode::deep_multilevel) {
            ASSERT(context.partition.objective != Objective::steiner_tree);
            partitioned_hypergraph = DeepMultilevel<TypeTraits>::partition(hypergraph, context);
        } else {
            throw InvalidParameterException("Invalid partitioning mode!");
        }

        Individual individual(partitioned_hypergraph, context);
        return population.addStartingIndividual(individual, context);
    } 

    template<typename TypeTraits>
    EvoDecision EvoPartitioner<TypeTraits>::decideNextMove(const Context& context) {
        if (utils::Randomize::instance().getRandomFloat(0, 1, THREAD_ID) < context.evolutionary.mutation_chance) {
            return EvoDecision::mutation;
        }
        return EvoDecision::combine;
    }

    template<typename TypeTraits>
    void EvoPartitioner<TypeTraits>::performCombine(EvoPartitioner<TypeTraits>::Hypergraph& hypergraph, const Context& context, TargetGraph* target_graph, Population& population) {
        //const auto& parents = population.tournamentSelect();
        //std::vector<PartitionID> part1 = parents.first.get().partition();
        //std::vector<PartitionID> part2 = parents.second.get().partition()
        const size_t position1 = population.randomIndividual();
        const size_t position2 = population.randomIndividualExcept(position1);
        std::vector<PartitionID> part1 = population.individualAt(position1).partition();
        std::vector<PartitionID> part2 = population.individualAt(position2).partition();
        vec<PartitionID> comms(hypergraph.initialNumNodes());
        PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);
        for ( const HypernodeID& hn : hypergraph.nodes() ) {
            /*if (part1[hn] == part2[hn]) {
                partitioned_hypergraph.setOnlyNodePart(hn, part1[hn]);
            } else {
                partitioned_hypergraph.setOnlyNodePart(hn, context.partition.k + hn);
            }*/
            //partitioned_hypergraph.setOnlyNodePart(hn, part1[hn] + context.partition.k * part2[hn]);
            comms[hn] = part1[hn] + context.partition.k * part2[hn];
        }
        hypergraph.setCommunityIDs(std::move(comms));
        //partitioned_hypergraph.initializePartition();
        if (context.partition.mode == Mode::direct) {
            //Multilevel<TypeTraits>::partitionVCycle(hypergraph, partitioned_hypergraph, context, target_graph);
            partitioned_hypergraph = Multilevel<TypeTraits>::partition(hypergraph, context, target_graph);
        } else {
            throw InvalidParameterException("Invalid partitioning mode!");
        }

        Individual individual(partitioned_hypergraph, context);    
        population.insert(std::move(individual), context);
    }

    // REMEMBER: When switching barck to calling fpartitionVCycleorce inserts, add the position i.e:
    // _population.forceInsertSaveBest(mutate::vCycleWithNewInitialPartitioning(hg, _population,_population.individualAt(mutation_position), context),  mutation_position);
    template<typename TypeTraits>
    void EvoPartitioner<TypeTraits>::performMutation(EvoPartitioner<TypeTraits>::Hypergraph& hypergraph, const Context& context, TargetGraph* target_graph, Population& population) {
        const size_t mutation_position = population.randomIndividual();
        std::vector<PartitionID> cur = population.individualAt(mutation_position).partition();
        PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);
        vec<PartitionID> comms(hypergraph.initialNumNodes());
        for ( const HypernodeID& hn : hypergraph.nodes() ) {
            //partitioned_hypergraph.setOnlyNodePart(hn, cur[hn]);
            comms[hn] = cur[hn];
        }
        hypergraph.setCommunityIDs(std::move(comms));
        //partitioned_hypergraph.initializePartition();
        if (context.partition.mode == Mode::direct) {
            //Multilevel<TypeTraits>::partitionVCycle(hypergraph, partitioned_hypergraph, context, target_graph);
            partitioned_hypergraph = Multilevel<TypeTraits>::partition(hypergraph, context, target_graph);
        } else {
            throw InvalidParameterException("Invalid partitioning mode!");
        }

        Individual individual(partitioned_hypergraph, context);    
        population.insert(std::move(individual), context);
    }

    template<typename TypeTraits>
    void EvoPartitioner<TypeTraits>::performEvolution(EvoPartitioner<TypeTraits>::Hypergraph& hg, Context& context, TargetGraph* target_graph, Population& population) { 
        int timelimit = context.partition.time_limit;
        utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
        int mutations = 0;
        int combinations = 0;
        while (timer.get("evolutionary") <= timelimit) {
            //LOG << timer.get("evolutionary") << " of " << timelimit << "\n";
            ++context.evolutionary.iteration;

            timer.start_timer("evolutionary", "Evolutionary");

            EvoDecision decision = decideNextMove(context);
            switch (decision) {
                case EvoDecision::mutation:
                    performMutation(hg, context, target_graph, population);
                    mutations++;
                    break;
                case EvoDecision::combine:
                    performCombine(hg, context, target_graph, population);
                    combinations++;
                    break;
                default:
                    LOG << "Error in evo_partitioner.cpp: Non-covered case in decision making";
                    std::exit(EXIT_FAILURE);
            }
            timer.stop_timer("evolutionary");
        }
        hg.reset();
        LOG << "Performed " << context.evolutionary.iteration << " Evolutionary Iterations" << "\n";
        LOG << "    " << mutations << " Mutations" << "\n";
        LOG << "    " << combinations << " Combinations" << "\n";
    }

}  // namespace mt_kahypar
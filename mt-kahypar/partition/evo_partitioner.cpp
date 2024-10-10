#include "evo_partitioner.h"
#include "partitioner.cpp"
#include <mutex>
#include <limits>

namespace mt_kahypar {

    template<typename TypeTraits>
    typename EvoPartitioner<TypeTraits>::PartitionedHypergraph EvoPartitioner<TypeTraits>::partition(
        const Hypergraph& hypergraph, Context& context, TargetGraph* target_graph) {
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
        /*DegreeZeroHypernodeRemover<TypeTraits> degree_zero_hn_remover(context);
        LargeHyperedgeRemover<TypeTraits> large_he_remover(context);
        preprocess(hypergraph, context, target_graph);
        sanitize(hypergraph, context, degree_zero_hn_remover, large_he_remover);*/
        timer.stop_timer("preprocessing");

        Population best_population;
        std::string history;
        for (size_t i = 0; i < context.evolutionary.repetitions; ++i) {
            Population population;
            context.evolutionary.iteration = 0;
            history += generateInitialPopulation(hypergraph, context, target_graph, population);

            history += performEvolution(hypergraph, context, target_graph, population);
            history += "\nRUN COMPLETED\n";
            // if (best_population.size() == 0 || best_population.bestFitness() > population.bestFitness()) {
            //     best_population = std::move(population);
            // }
            best_population.merge(population, context.evolutionary.output_size);
        }

        if (context.evolutionary.frequency_file != "") {
            if constexpr (Hypergraph::is_graph) {
                std::vector<uint32_t> frequencies = best_population.getEdgeFrequencies(hypergraph, false);
                std::ofstream out_stream(context.evolutionary.frequency_file.c_str());
                out_stream << "# max=" << context.evolutionary.output_size << std::endl;
                out_stream << "id_high_degree,id_low_degree,frequency" << std::endl;
                for (HyperedgeID edge: hypergraph.edges()) {
                    HypernodeID u = hypergraph.edgeSource(edge);
                    HypernodeID v = hypergraph.edgeTarget(edge);
                    HyperedgeID deg_u = hypergraph.nodeDegree(u);
                    HyperedgeID deg_v = hypergraph.nodeDegree(v);
                    if (deg_u > deg_v || (deg_u == deg_v && u >= v)) {
                        out_stream << u << "," << v << "," << frequencies[edge] << std::endl;
                    }
                }
                out_stream.close();
            } else {
                WARNING("Outputting edge frequencies for hypergraphs is not implemented!");
            }
        }

        LOG << "\nFinal Best:   " << best_population.individualAt(best_population.best()).fitness();
        LOG << "Final Worst:  " << best_population.individualAt(best_population.worst()).fitness();
        LOG << "Final Size:   " << best_population.size() << "\n";

        if (context.evolutionary.history_file != "") {
            std::ofstream out_stream(context.evolutionary.history_file.c_str());
            out_stream << history;
            out_stream.close();
        }

        PartitionedHypergraph partitioned_hypergraph(context.partition.k, const_cast<Hypergraph&>(hypergraph));
        
        std::vector<PartitionID> best = best_population.individualAt(best_population.best()).partition();
        partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
            partitioned_hypergraph.setOnlyNodePart(hn, best[hn]);
        });
        partitioned_hypergraph.initializePartition();

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
        /*large_he_remover.restoreLargeHyperedges(partitioned_hypergraph);
        degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hypergraph);
        forceFixedVertexAssignment(partitioned_hypergraph, context);*/
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
    std::string EvoPartitioner<TypeTraits>::generateInitialPopulation(const Hypergraph& hg, Context& context, TargetGraph* target_graph, Population& population) {
        context.partition.verbose_output = false;
        int timelimit = context.partition.time_limit;
        utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
        auto start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
        auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
        auto time_elapsed = now - start;
        auto duration = std::chrono::seconds(timelimit);
        std::string history = "Starttime: " + std::to_string(start.count()) + "\n";
        // INITIAL POPULATION
        if (context.evolutionary.dynamic_population_size) {
            timer.start_timer("evolutionary", "Evolutionary");
            auto fitness = generateIndividual(hg, context, target_graph, population).fitness();
            now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
            history += "" + std::to_string(now.count()) + ", Initial, " + std::to_string(fitness) + "\n";
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
        HyperedgeWeight best = std::numeric_limits<HyperedgeWeight>::max();
        int iteration = 0;
        while (population.size() < context.evolutionary.population_size && 
            time_elapsed <= duration) {
            ++context.evolutionary.iteration;
            timer.start_timer("evolutionary", "Evolutionary");
            auto cur = generateIndividual(hg, context, target_graph, population).fitness();
            timer.stop_timer("evolutionary");
            now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
            if (iteration == 0 || (cur < best)) {
                best = cur;
                history += "" + std::to_string(now.count()) + ", Repetition, " + std::to_string(cur) + "\n";
            }
            iteration++;
            time_elapsed = now - start;
        }
        context.evolutionary.time_elapsed = time_elapsed;
        context.partition.verbose_output = true;
        return history;
    }

    template<typename TypeTraits>
    const Individual & EvoPartitioner<TypeTraits>::generateIndividual(const Hypergraph& input_hg, Context& context, TargetGraph* target_graph, Population& population) {
        Hypergraph hypergraph = input_hg.copy(parallel_tag_t{});
        DegreeZeroHypernodeRemover<TypeTraits> degree_zero_hn_remover(context);
        LargeHyperedgeRemover<TypeTraits> large_he_remover(context);
        preprocess(hypergraph, context, target_graph);
        sanitize(hypergraph, context, degree_zero_hn_remover, large_he_remover);
        
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

        large_he_remover.restoreLargeHyperedges(partitioned_hypergraph);
        degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hypergraph);
        forceFixedVertexAssignment(partitioned_hypergraph, context);

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
    EvoMutateStrategy EvoPartitioner<TypeTraits>::decideNextMutation(const Context& context) {
        if (utils::Randomize::instance().flipCoin(THREAD_ID)) {
            return EvoMutateStrategy::vcycle;
        }
        return EvoMutateStrategy::new_initial_partitioning_vcycle;
    }


    template<typename TypeTraits>
    vec<PartitionID> EvoPartitioner<TypeTraits>::combinePartitions(const Context& context, Population& population, const std::vector<size_t>& ids) {
        vec<PartitionID> combined(population.individualAt(0).partition().size());

        std::unordered_map<std::string, int> tuple_to_block;
        int current_community = 0;

        for (size_t vertex = 0; vertex < combined.size(); vertex++) {
            std::string partition_tuple;
            for (auto id : ids) {
                partition_tuple += std::to_string(population.individualAt(id).partition()[vertex]) + ",";
            }

            if (tuple_to_block.find(partition_tuple) == tuple_to_block.end()) {
                tuple_to_block[partition_tuple] = current_community++;
            }

            combined[vertex] = tuple_to_block[partition_tuple];
        }

        return combined;
    }

    template<typename TypeTraits>
    std::string EvoPartitioner<TypeTraits>::performCombine(const Hypergraph& input_hg, const Context& context, TargetGraph* target_graph, Population& population) {
        std::vector<size_t> parents;
        size_t best = population.randomIndividual();
        parents.push_back(best);
        for (int x = 1; x < context.evolutionary.kway_combine; x ++) {
            size_t new_parent = population.randomIndividual();
            parents.push_back(new_parent);
            if (population.individualAt(new_parent).fitness() <= population.individualAt(best).fitness()) {
                best = new_parent;
            }
        }
        std::unordered_map<PartitionID, int> comm_to_block;
        vec<PartitionID> comms = combinePartitions(context, population, parents);
        vec<PartitionID> part(input_hg.initialNumNodes());

        Hypergraph hypergraph = input_hg.copy(parallel_tag_t{});
        PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);
        for ( const HypernodeID& hn : hypergraph.nodes() ) { 
            partitioned_hypergraph.setOnlyNodePart(hn, population.individualAt(best).partition()[hn]);
            if (comm_to_block.find(comms[hn]) == comm_to_block.end()) {
                comm_to_block[comms[hn]] = population.individualAt(best).partition()[hn];
            }
        }
        partitioned_hypergraph.initializePartition();
        hypergraph.setCommunityIDs(std::move(comms));
        if (context.partition.mode == Mode::direct) {
            Multilevel<TypeTraits>::evolutionPartitionVCycle(hypergraph, partitioned_hypergraph, context, comm_to_block, target_graph);
        } else {
            throw InvalidParameterException("Invalid partitioning mode!");
        }

        Individual individual(partitioned_hypergraph, context);
        //LOG << "Combined Individuals with fitness: " << population.individualAt(position1).fitness() << "," << population.individualAt(position2).fitness() << " into " << individual.fitness();
        std::string ret = "";
        if (individual.fitness() < population.bestFitness()) {
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
            ret = "" + std::to_string(time.count()) + ", Combine, " + std::to_string(individual.fitness()) + "\n";
        }
        std::lock_guard<std::mutex> lock = population.getLock();
        population.insert(std::move(individual), context);
        return ret;
    }

    template<typename TypeTraits>
    std::string EvoPartitioner<TypeTraits>::performMutation(const Hypergraph& input_hg, const Context& context, TargetGraph* target_graph, Population& population) {
        Hypergraph hypergraph = input_hg.copy(parallel_tag_t{});
        const size_t mutation_position = population.randomIndividual();
        std::vector<PartitionID> cur = population.individualAt(mutation_position).partition();
        PartitionedHypergraph partitioned_hypergraph(context.partition.k, hypergraph);
        EvoMutateStrategy mutation = decideNextMutation(context);
        if (mutation == EvoMutateStrategy::vcycle) {
            //LOG << "Mutating old: ";
            vec<PartitionID> comms(hypergraph.initialNumNodes());
            std::unordered_map<PartitionID, int> comm_to_block;
            for ( const HypernodeID& hn : hypergraph.nodes() ) {
                partitioned_hypergraph.setOnlyNodePart(hn, cur[hn]);
                comms[hn] = cur[hn];
            }
            for (PartitionID i = 0; i < context.partition.k; i++) {
                comm_to_block[i] = i;
            }
            partitioned_hypergraph.initializePartition();
            hypergraph.setCommunityIDs(std::move(comms));
            if (context.partition.mode == Mode::direct) {
                Multilevel<TypeTraits>::evolutionPartitionVCycle(hypergraph, partitioned_hypergraph, context, comm_to_block, target_graph);
            } else {
                throw InvalidParameterException("Invalid partitioning mode!");
            }
        } else if (mutation == EvoMutateStrategy::new_initial_partitioning_vcycle) {
            //LOG << "Mutating new: ";
            vec<PartitionID> comms(hypergraph.initialNumNodes());
            for ( const HypernodeID& hn : hypergraph.nodes() ) {
                comms[hn] = cur[hn];
            }
            hypergraph.setCommunityIDs(std::move(comms));
            if (context.partition.mode == Mode::direct) {
                partitioned_hypergraph = Multilevel<TypeTraits>::partition(hypergraph, context, target_graph);
            } else {
                throw InvalidParameterException("Invalid partitioning mode!");
            }
        }
        Individual individual(partitioned_hypergraph, context);
        //LOG << "Mutated Individual with fitness: " << population.individualAt(mutation_position).fitness() << " to " << individual.fitness();
        std::string ret = "";
        if (individual.fitness() < population.bestFitness()) {
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
            if (mutation == EvoMutateStrategy::vcycle) {
                ret = "" + std::to_string(time.count()) + ", MutateOld, " + std::to_string(individual.fitness()) + "\n";
            } else {
                ret = "" + std::to_string(time.count()) + ", MutateNew, " + std::to_string(individual.fitness()) + "\n";
            }
        }
        std::lock_guard<std::mutex> lock = population.getLock();
        population.insert(std::move(individual), context);
        return ret;
    }

    inline void disableTimerAndStatsEvo(const Context& context) {
        if ( context.type == ContextType::main && context.partition.mode == Mode::direct ) {
            utils::Utilities& utils = utils::Utilities::instance();
            parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
            utils.getTimer(context.utility_id).disable();
            utils.getStats(context.utility_id).disable();
        }
    }

    inline void enableTimerAndStatsEvo(const Context& context) {
        if ( context.type == ContextType::main && context.partition.mode == Mode::direct ) {
            utils::Utilities& utils = utils::Utilities::instance();
            parallel::MemoryPool::instance().activate_unused_memory_allocations();
            utils.getTimer(context.utility_id).enable();
            utils.getStats(context.utility_id).enable();
        }
    }

    template<typename TypeTraits>
    std::string EvoPartitioner<TypeTraits>::performEvolution(const Hypergraph& hg, Context& context, TargetGraph* target_graph, Population& population) {
        context.partition.verbose_output = false;
        int timelimit = context.partition.time_limit;
        utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id);
        int mutations = 0;
        int combinations = 0;
        auto time_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
        std::string history = "";
        std::mutex _history_mutex;
        auto duration = std::chrono::seconds(timelimit) - context.evolutionary.time_elapsed;
        std::atomic<bool> stop_flag(false);
        timer.start_timer("evolutionary", "Evolutionary");
        while(!stop_flag) {
            auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch());
            if (now - time_start >= duration) {
                stop_flag = true;
                break;
            }
            Context evo_context(context);
            evo_context.type = ContextType::main;
            evo_context.utility_id = utils::Utilities::instance().registerNewUtilityObjects();
            EvoDecision decision = decideNextMove(context);
            EvoPartitioner<TypeTraits>::Hypergraph hg_copy = hg.copy();
            switch (decision) {
                case EvoDecision::mutation: 
                    {
                        std::string h = performMutation(hg_copy, evo_context, target_graph, population);
                        std::lock_guard<std::mutex> lock(_history_mutex);
                        history += h;
                        mutations++;
                        ++context.evolutionary.iteration;
                        break;
                    }
                case EvoDecision::combine: 
                    {
                        std::string h = performCombine(hg_copy, evo_context, target_graph, population);
                        combinations++;
                        std::lock_guard<std::mutex> lock(_history_mutex);
                        history += h;
                        ++context.evolutionary.iteration;
                        break;
                    }
                default:
                    LOG << "Error in evo_partitioner.cpp: Non-covered case in decision making";
                    std::exit(EXIT_FAILURE);
            }
        }
        timer.stop_timer("evolutionary");

        context.partition.verbose_output = true;
        LOG << "\nPerformed " << context.evolutionary.iteration << " Evolutionary Iterations";
        LOG << "    " << (context.evolutionary.iteration - mutations - combinations)
            << " Initial Population members";
        LOG << "    " << mutations << " Mutations";
        LOG << "    " << combinations << " Combinations";
        LOG << "Best:   " << population.individualAt(population.best()).fitness();
        LOG << "Half:   " << population.listOfBest((population.size() + 1) / 2).back().get().fitness();
        LOG << "Worst:  " << population.individualAt(population.worst()).fitness() << "\n";

        return history;
    }

}  // namespace mt_kahypar
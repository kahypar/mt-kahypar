#include "mt-kahypar/partition/multicandidate.h"

#include <memory>

#include "include/libmtkahypartypes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/coarsening/multilevel/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/hypergraph_statistics.h"
#include "mt-kahypar/utils/utilities.h"
#include "tbb/parallel_for.h"


#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"
#include "mt-kahypar/partition/deep_multilevel.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/partition/coarsening/nlevel/nlevel_uncoarsener.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {

namespace {

void disableTimerAndStats(const Context& context) {
  if (context.type == ContextType::main &&
      context.partition.mode == Mode::direct) {
    utils::Utilities& utils = utils::Utilities::instance();
    parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
    utils.getTimer(context.utility_id).disable();
    utils.getStats(context.utility_id).disable();
  }
}

void enableTimerAndStats(const Context& context) {
  if (context.type == ContextType::main &&
      context.partition.mode == Mode::direct) {
    utils::Utilities& utils = utils::Utilities::instance();
    parallel::MemoryPool::instance().activate_unused_memory_allocations();
    utils.getTimer(context.utility_id).enable();
    utils.getStats(context.utility_id).enable();
  }
}

  template <typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph multicandidate_partitioning(
      typename TypeTraits::Hypergraph& hypergraph,
      const Context& context, 
      const TargetGraph* target_graph) {
    using Hypergraph = typename TypeTraits::Hypergraph;
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    // ################## COARSENING ##################
    mt_kahypar::io::printCoarseningBanner(context);
  
    const bool nlevel = context.isNLevelPartitioning();
    UncoarseningData<TypeTraits> uncoarseningData(nlevel, hypergraph, context);
    utils::Timer& timer = utils::Utilities::instance().getTimer(context.utility_id); 
    timer.start_timer("coarsening", "Coarsening");
    {
      std::unique_ptr<ICoarsener> coarsener =
          CoarsenerFactory::getInstance().createObject(
              context.coarsening.algorithm, utils::hypergraph_cast(hypergraph),
              context, uncoarsening::to_pointer(uncoarseningData));
      coarsener->coarsen();
      if (context.partition.verbose_output) {
        mt_kahypar_hypergraph_t coarsestHypergraph = coarsener->coarsestHypergraph();
        mt_kahypar::io::printHypergraphInfo(
          utils::cast<Hypergraph>(coarsestHypergraph), context,
          "Coarsened Hypergraph", context.partition.show_memory_consumption);
      }
    }
    timer.stop_timer("coarsening");
  
    // ################## INITIAL PARTITIONING ##################

    PartitionedHypergraph& partitioned_hg = uncoarseningData.coarsestPartitionedHypergraph();
    vec<vec<PartitionID>> partition_pool;
    auto uncoarsener = std::make_unique<MultilevelUncoarsener<TypeTraits>>(
        hypergraph, context, uncoarseningData, target_graph);

    auto ip = [&]() {
      timer.start_timer("initial_partitioning", "Initial Partitioning");
      partitioned_hg.resetData();
      uncoarsener->updateMetrics();
      //GainCachePtr::resetGainCache(uncoarsener->getGainCache());
      DegreeZeroHypernodeRemover<TypeTraits> degree_zero_hn_remover(context);
      if (context.initial_partitioning.remove_degree_zero_hns_before_ip) {
        degree_zero_hn_remover.removeDegreeZeroHypernodes(
            partitioned_hg.hypergraph());
      }

      Context ip_context(context);
      ip_context.type = ContextType::initial_partitioning;
      ip_context.refinement = context.initial_partitioning.refinement;
      disableTimerAndStats(context);
      if (context.initial_partitioning.mode == Mode::direct) {
        // The pool initial partitioner consist of several flat bipartitioning
        // techniques. This case runs as a base case (k = 2) within recursive
        // bipartitioning or the deep multilevel scheme.
        ip_context.partition.verbose_output = false;
        Pool<TypeTraits>::bipartition(partitioned_hg, ip_context);
      } else if (context.initial_partitioning.mode ==
                 Mode::recursive_bipartitioning) {
        RecursiveBipartitioning<TypeTraits>::partition(
            partitioned_hg, ip_context, target_graph);
      } else if (context.initial_partitioning.mode == Mode::deep_multilevel) {
        ASSERT(ip_context.partition.objective != Objective::steiner_tree);
        ip_context.partition.verbose_output = false;
        DeepMultilevel<TypeTraits>::partition(partitioned_hg, ip_context);
      } else {
        throw InvalidParameterException(
            "Undefined initial partitioning algorithm");
      }
      enableTimerAndStats(context);
      degree_zero_hn_remover.restoreDegreeZeroHypernodes(partitioned_hg);
      //io::printPartitioningResults(partitioned_hg, context,
      //                             "Initial Partitioning Results:");
      uncoarsener->updateMetrics();
      //GainCachePtr::resetGainCache(uncoarsener->getGainCache());
      uncoarsener->refine();
      partition_pool.push_back({});
      partition_pool.back().resize(partitioned_hg.initialNumNodes());
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        partition_pool.back()[hn] = partitioned_hg.partID(hn);
      });
      timer.stop_timer("initial_partitioning");
    };

    // ################## UNCOARSENING ##################

    // Functions

    auto replacePartition = [&partitioned_hg, &uncoarsener, &context](const vec<PartitionID>& partition) {
      partitioned_hg.resetData();
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        const PartitionID part_id = partition[hn];
        partitioned_hg.setOnlyNodePart(hn, part_id);
      });
      partitioned_hg.initializePartition();
      uncoarsener->updateMetrics();
      GainCachePtr::resetGainCache(uncoarsener->getGainCache());
    };

    // Maybe replace with comparison of PartitioningResult
    auto isBetterThan = [&partitioned_hg, replacePartition](
                            const vec<PartitionID>& l,
                            const vec<PartitionID>& r) {
      vec<PartitionID> tmp_ids;
      tmp_ids.resize(partitioned_hg.initialNumNodes());
      partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
        tmp_ids[hn] = partitioned_hg.partID(hn);
      });

      replacePartition(l);
      auto l_quality = metrics::quality(partitioned_hg, Objective::km1);

      replacePartition(r);
      auto r_quality = metrics::quality(partitioned_hg, Objective::km1);

      replacePartition(tmp_ids);
      return l_quality < r_quality;
    };
    //

    io::printLocalSearchBanner(context);
    //timer.start_timer("refinement", "Refinement");
    uncoarsener->initialize();
    uncoarsener->stepNextLevel();
    int level = 3; //TODO: make this a parameter
    while(!uncoarsener->isTopLevel()) {
      if (level > 0) {
        for (int i = 0; i < level * 5; i++) {
          ip();
        }

        std::cout << "level " << level << std::endl;
        std::cout << "Number of partitions: " << partition_pool.size() << std::endl;

        // Reset phg to contain the first partition of the partition pool
        replacePartition(partition_pool[0]);
        uncoarsener->projectToNextLevelAndRefine();

        partition_pool[0].resize(partitioned_hg.initialNumNodes());
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
          partition_pool[0][hn] = partitioned_hg.partID(hn);
        });
        // Project the partitions in the partition pool to the next level
        for(auto it = std::next(partition_pool.begin()); it != partition_pool.end(); it++) {
          vec<PartitionID> tmp_partition;
          tmp_partition.resize(partitioned_hg.initialNumNodes());
          std::copy(it->begin(), it->end(), tmp_partition.begin());

          it->resize(partitioned_hg.initialNumNodes());
          partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
            const HypernodeID coarse_hn = (uncoarseningData.hierarchy)[uncoarsener->currentLevel() + 1].mapToContractedHypergraph(hn);
            const PartitionID part_id = tmp_partition[coarse_hn];
            (*it)[hn] = part_id;
          });
        }
        // Refine all the partitions
        for (auto it = std::next(partition_pool.begin()); it != partition_pool.end(); it++) {
          replacePartition(*it);
          uncoarsener->refine();
          // Write back the refined partition
          it->resize(partitioned_hg.initialNumNodes());
          partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {

            (*it)[hn] = partitioned_hg.partID(hn);
          });
        }

        // Do local tournament
        if(level > 1) {
          if(partition_pool.size() <= 2) {
            --level;
            continue;
          }
          // Sort the array of partitioned hypergraphs by quality
          std::sort(partition_pool.begin(), partition_pool.end(), isBetterThan);
          // first step: mark the winners of each match
          // mark the first 1 or 2 elements as winners depending on whether 
          // the number of partitions is even or odd
          std::vector<bool> winners(partition_pool.size(), true);
          int offset = partition_pool.size() % 2 ? 1 : 2;
          utils::Randomize::instance().parallelShuffleVector(partition_pool, offset, partition_pool.size());
          for(size_t i = offset; i < partition_pool.size(); i += 2) {
            if(isBetterThan(partition_pool[i], partition_pool[i + 1])) {
              winners[i + 1] = false;
            } else {
              winners[i] = false;
            }
          }
          //erase the losers
          partition_pool.erase(std::remove_if(partition_pool.begin(), partition_pool.end(),
                                              [&](const vec<PartitionID>& p) {
                                                return !winners[&p - &partition_pool[0]];
                                              }),
                              partition_pool.end());
        } else {
          int round = 0;
          // Do final tournament and choose winner and set the partitioned hypergraph to the winner
          while(partition_pool.size() > 1 && !uncoarsener->isTopLevel()) {
            std::cout << "Round " << round++ << std::endl;
            std::cout << "Number of partitions: " << partition_pool.size() << std::endl;
            // Print qualities of the partitions
            for(auto it = partition_pool.begin(); it != partition_pool.end(); it++) {
              replacePartition(*it);
              auto quality = metrics::quality(partitioned_hg, Objective::km1);
              std::cout << "Quality: " << quality << std::endl;
            }
            int offset = partition_pool.size() % 2 ? 1 : 0;
            utils::Randomize::instance().parallelShuffleVector(
                partition_pool, offset, partition_pool.size());

            std::vector<bool> winners(partition_pool.size(), true);
            for (size_t i = offset; i < partition_pool.size(); i += 2) {
              if (isBetterThan(partition_pool[i], partition_pool[i + 1])) {
                winners[i + 1] = false;
              } else {
                winners[i] = false;;
              }
            }
            // erase the losers
            partition_pool.erase(
                std::remove_if(partition_pool.begin(), partition_pool.end(),
                               [&](const vec<PartitionID>& p) {
                                 return !winners[&p - &partition_pool[0]];
                               }),
                partition_pool.end());
            // Project the partitions in the partition pool to the next level and refine
            replacePartition(partition_pool[0]);
            uncoarsener->projectToNextLevelAndRefine();
            partition_pool[0].resize(partitioned_hg.initialNumNodes());
            partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
              partition_pool[0][hn] = partitioned_hg.partID(hn);
            });
            for(auto it = std::next(partition_pool.begin()); it != partition_pool.end(); it++) {
              vec<PartitionID> tmp_partition;
              tmp_partition.resize(partitioned_hg.initialNumNodes());
              std::copy(it->begin(), it->end(), tmp_partition.begin());

              it->resize(partitioned_hg.initialNumNodes());
              partitioned_hg.doParallelForAllNodes([&](const HypernodeID hn) {
                const HypernodeID coarse_hn = (uncoarseningData.hierarchy)[uncoarsener->currentLevel() + 1].mapToContractedHypergraph(hn);
                const PartitionID part_id = tmp_partition[coarse_hn];
                (*it)[hn] = part_id;
              });
            }
            for(auto it = partition_pool.begin(); it != partition_pool.end(); it++) {
              replacePartition(*it);
              uncoarsener->refine();
              // Write back the refined partition
              ds::Array<PartitionID> tmp_ids;
              tmp_ids.resize(partitioned_hg.initialNumNodes());
              partitioned_hg.extractPartIDs(tmp_ids);
              std::copy(tmp_ids.begin(), tmp_ids.end(), (*it).begin());
            }
          }
          std::sort(partition_pool.begin(), partition_pool.end(), isBetterThan);
          replacePartition(partition_pool[0]);
          std::cout << "Final partition quality: " << metrics::quality(partitioned_hg, Objective::km1) << std::endl;
        }
        --level;
      } else {
        uncoarsener->projectToNextLevelAndRefine();
      }
    }
    if(partition_pool.size() == 0) {
      ip();
    } else if (partition_pool.size() > 1){
      std::sort(partition_pool.begin(), partition_pool.end(), isBetterThan);
      replacePartition(partition_pool[0]);
    }
    uncoarsener->rebalancing();
    io::printPartitioningResults(partitioned_hg, context, "Local Search Results:");
    return uncoarsener->movePartitionedHypergraph();
  }
} // namespace

template<typename TypeTraits>
typename Multicandidate<TypeTraits>::PartitionedHypergraph
Multicandidate<TypeTraits>::partition(Hypergraph& hypergraph,
                                      const Context& context,
                                      const TargetGraph* target_graph) {
  PartitionedHypergraph partitioned_hg = 
    multicandidate_partitioning<TypeTraits>(hypergraph, context, target_graph);
  return partitioned_hg;
}

template<typename TypeTraits>
void Multicandidate<TypeTraits>::partition(
    PartitionedHypergraph& partitioned_hg, const Context& context,
    const TargetGraph* target_graph) {
  PartitionedHypergraph tmp_phg =
      partition(partitioned_hg.hypergraph(), context, target_graph);
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(Multicandidate)

} // namespace mt_kahypar
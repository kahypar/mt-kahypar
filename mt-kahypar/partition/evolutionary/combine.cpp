#include "mt-kahypar/partition/evolutionary/combine.h"


namespace mt_kahypar::combine {
  vec<PartitionID> combinePartitions(const std::vector<std::vector<PartitionID>> &parent_partitions) {
    // NM: ugly implementation, alternative suggestion: use something like MSD radix sort, i.e., split in k
    // buckets and recurse on each non-empty bucket, with each recursion level corresponding to a different parent
    // partition. Return first non-used index so it is possible to still guarantee contiguous community indices
    // (low prio)

    vec<PartitionID> combined(parent_partitions[0].size());
    std::unordered_map<std::string, int> tuple_to_block;
    int current_community = 0;

    for (int vertex = 0; vertex < combined.size(); vertex++) {
      std::string partition_tuple;
      for (size_t i = 0; i < parent_partitions.size(); ++i) {
        partition_tuple += std::to_string(parent_partitions[i][vertex]) + ",";
      }

      if (tuple_to_block.find(partition_tuple) == tuple_to_block.end()) {
        tuple_to_block[partition_tuple] = current_community++;
      }

      combined[vertex] = tuple_to_block[partition_tuple];
    }

    return combined;
  }


  Context modifyContext(const Context& context, ContextModifierParameters params) {
    Context modifiedContext(context);
    modifiedContext.partition.k = params.k;
    modifiedContext.partition.epsilon = params.epsilon;
    modifiedContext.partition.mode = params.recursive_bipartitioning ? Mode::recursive_bipartitioning : modifiedContext.partition.mode;
    return modifiedContext;
  }



}


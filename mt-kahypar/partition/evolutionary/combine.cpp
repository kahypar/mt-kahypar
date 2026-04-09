#include "mt-kahypar/partition/evolutionary/combine.h"


namespace mt_kahypar::combine {
  vec<PartitionID> combinePartitions(mt_kahypar::Population& population, const std::vector<size_t>& ids) {
    // aquire lock --- possibly unnecessary
    std::vector<std::vector<PartitionID>> parent_partitions;
    for (auto id : ids) {
      parent_partitions.push_back(population.partitionCopySafe(id)); // FIXED
    }

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

}


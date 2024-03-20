#pragma once

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

// Forward Declaration
class TargetGraph;

template <typename TypeTraits>
class Multicandidate {
  
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  static PartitionedHypergraph partition(Hypergraph& hypergraph,
                                         const Context& context,
                                         const TargetGraph* target_graph = nullptr);
                                         
  static void partition(PartitionedHypergraph& hypergraph,
                        const Context& context,
                        const TargetGraph* target_graph = nullptr);
};

}  // namespace mt_kahypar

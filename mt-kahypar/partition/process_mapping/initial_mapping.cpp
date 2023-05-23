/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/process_mapping/initial_mapping.h"

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/process_mapping/process_graph.h"
#include "mt-kahypar/partition/process_mapping/dual_bipartitioning.h"
#include "mt-kahypar/partition/process_mapping/greedy_mapping.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

namespace {

using HyperedgeVector = vec<vec<HypernodeID>>;

template<typename TargetType, typename PartitionedHypergraph>
typename TargetType::PartitionedHypergraph convert(const PartitionedHypergraph& phg) {
  using Hypergraph = typename TargetType::Hypergraph;
  using TargetPartitionedHypergraph = typename TargetType::PartitionedHypergraph;
  using Factory = typename Hypergraph::Factory;

  const HypernodeID num_hypernodes = phg.initialNumNodes();
  const HyperedgeID num_hyperedges = phg.initialNumEdges();
  HyperedgeVector edge_vector;
  vec<HyperedgeWeight> hyperedge_weight;
  vec<HypernodeWeight> hypernode_weight;

  // Allocate data structure
  tbb::parallel_invoke([&] {
    edge_vector.assign(num_hyperedges, vec<HypernodeID>());
  }, [&] {
    hyperedge_weight.assign(num_hyperedges, 0);
  }, [&] {
    hypernode_weight.assign(num_hypernodes, 0);
  });

  // Write hypergraph into temporary data structure
  tbb::parallel_invoke([&] {
    phg.doParallelForAllEdges([&](const HyperedgeID& he) {
      hyperedge_weight[he] = phg.edgeWeight(he);
      for ( const HypernodeID& pin : phg.pins(he) ) {
        edge_vector[he].push_back(pin);
      }
    });
  }, [&] {
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      hypernode_weight[hn] = phg.nodeWeight(hn);
    });
  });

  // Construct new hypergraph and apply partition
  Hypergraph converted_hg = Factory::construct(num_hypernodes, num_hyperedges,
    edge_vector, hyperedge_weight.data(), hypernode_weight.data());
  TargetPartitionedHypergraph converted_phg(phg.k(), converted_hg, parallel_tag_t { });
  phg.doParallelForAllNodes([&](const HypernodeID& hn) {
    converted_phg.setOnlyNodePart(hn, phg.partID(hn));
  });
  converted_phg.initializePartition();

  ASSERT(metrics::quality(phg, Objective::cut) == metrics::quality(converted_phg, Objective::cut));
  return converted_phg;
}

template<typename PartitionedHypergraph, typename TargetPartitionedHypergraph>
void applyPartition(PartitionedHypergraph& phg,
                    TargetPartitionedHypergraph& target_phg) {
  target_phg.doParallelForAllNodes([&](const HypernodeID& hn) {
    const PartitionID from = target_phg.partID(hn);
    const PartitionID to = phg.partID(hn);
    if ( from != to ) {
      target_phg.changeNodePart(hn, from, to);
    }
  });
}

template<typename PartitionedHypergraph>
void map_to_process_graph(PartitionedHypergraph& communication_hg,
                          const ProcessGraph& process_graph,
                          const Context& context) {
  using Hypergraph = typename PartitionedHypergraph::UnderlyingHypergraph;
  // We contract all blocks of the partition to create an one-to-one mapping problem
  vec<HypernodeID> mapping(communication_hg.initialNumNodes(), kInvalidHypernode);
  communication_hg.setProcessGraph(&process_graph);
  communication_hg.doParallelForAllNodes([&](const HypernodeID hn) {
    mapping[hn] = communication_hg.partID(hn);
  });
  // Here, we collapse each block of the communication hypergraph partition into
  // a single node. The contracted hypergraph has exactly k nodes. In the
  // contracted hypergraph node i corresponds to block i of the input
  // communication hypergraph.
  Hypergraph contracted_hg = communication_hg.hypergraph().contract(mapping);
  ASSERT(contracted_hg.initialNumNodes() == communication_hg.k());
  PartitionedHypergraph contracted_phg(communication_hg.k(), contracted_hg);
  for ( const HypernodeID& hn : contracted_phg.nodes() ) {
    contracted_phg.setOnlyNodePart(hn, hn);
  }
  contracted_phg.initializePartition();
  contracted_phg.setProcessGraph(&process_graph);

  const HyperedgeWeight objective_before = metrics::quality(contracted_phg, Objective::process_mapping);
  ASSERT(metrics::quality(communication_hg, Objective::process_mapping) == objective_before);

  // Solve one-to-one process mapping problem
  if ( context.process_mapping.strategy == ProcessMappingStrategy::dual_bipartitioning ) {
    DualBipartitioning<PartitionedHypergraph>::mapToProcessGraph(contracted_phg, process_graph, context);
  } else if ( context.process_mapping.strategy == ProcessMappingStrategy::greedy_mapping ) {
    GreedyMapping<PartitionedHypergraph>::mapToProcessGraph(contracted_phg, process_graph, context);
  }

  const HyperedgeWeight objective_after = metrics::quality(contracted_phg, Objective::process_mapping);
  if ( objective_after < objective_before ) {
    if ( context.partition.verbose_output ) {
      LOG << GREEN << "Initial process mapping algorithm has improved objective by"
          << (objective_before - objective_after)
          << "( Before =" << objective_before << ", After =" << objective_after << ")" << END;
    }
    // Initial mapping algorithm has improved solution quality
    // => apply improved mapping to input communication hypergraph
    communication_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID from = communication_hg.partID(hn);
      const PartitionID to = contracted_phg.partID(mapping[hn]);
      if ( from != to ) {
        communication_hg.changeNodePart(hn, from, to);
      }
    });
    ASSERT(metrics::quality(communication_hg, Objective::process_mapping) == objective_after);
  } else if ( context.partition.verbose_output && objective_before < objective_after ) {
    // Initial mapping algorithm has worsen solution quality
    // => use input partition of communication hypergraph
    LOG << RED << "Initial process mapping algorithm has worsen objective by"
      << (objective_after - objective_before)
      << "( Before =" << objective_before << ", After =" << objective_after << ")."
      << "Use mapping from initial partitiong!"<< END;
  }
}

}

template<typename TypeTraits>
void InitialMapping<TypeTraits>::mapToProcessGraph(PartitionedHypergraph& communication_hg,
                                                   const ProcessGraph& process_graph,
                                                   const Context& context) {
  if constexpr ( !PartitionedHypergraph::is_static_hypergraph ) {
    // The mapping algorithm collapses each block of the communication hypergraph partition into
    // a single node. Thereby, it uses the contract(...) function of the hypergraph data structure
    // which is only implemented for static graphs and hypergraphs (implemented for multilevel partitioing,
    // but not for n-level). In case the communication hypergraph uses an dynamic graph or hypergraph
    // data structure, we convert it to static data structure and then compute the initial mapping.
    if constexpr ( PartitionedHypergraph::is_graph ) {
      using TargetPartitionedGraph = typename StaticGraphTypeTraits::PartitionedHypergraph;
      TargetPartitionedGraph tmp_communication_hg = convert<StaticGraphTypeTraits>(communication_hg);
      map_to_process_graph(tmp_communication_hg, process_graph, context);
      applyPartition(tmp_communication_hg, communication_hg);
    } else {
      using TargetPartitionedHypergraph = typename StaticHypergraphTypeTraits::PartitionedHypergraph;
      TargetPartitionedHypergraph tmp_communication_hg = convert<StaticHypergraphTypeTraits>(communication_hg);
      map_to_process_graph(tmp_communication_hg, process_graph, context);
      applyPartition(tmp_communication_hg, communication_hg);
    }
  } else {
    map_to_process_graph(communication_hg, process_graph, context);
  }
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(InitialMapping)

}  // namespace kahypar

/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "tbb/parallel_for.h"

#include <string>
#include <vector>
#include <iostream>

#include "include/libmtkahypartypes.h"
#include "include/helper_functions.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/utils/cast.h"

namespace py = pybind11;
using namespace mt_kahypar;

namespace {
  void initialize_thread_pool(const size_t num_threads) {
    size_t P = num_threads;
    size_t num_available_cpus = mt_kahypar::HardwareTopology::instance().num_cpus();
    if ( num_available_cpus < num_threads ) {
      WARNING("There are currently only" << num_available_cpus << "cpus available."
        << "Setting number of threads from" << num_threads
        << "to" << num_available_cpus);
      P = num_available_cpus;
    }

    // Initialize TBB task arenas on numa nodes
    mt_kahypar::TBBInitializer::instance(P);
    // We set the membind policy to interleaved allocations in order to
    // distribute allocations evenly across NUMA nodes
    hwloc_cpuset_t cpuset = mt_kahypar::TBBInitializer::instance().used_cpuset();
    mt_kahypar::parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
    hwloc_bitmap_free(cpuset);
  }

  template<typename PartitionedHypergraph>
  double imbalance(const PartitionedHypergraph& partitioned_graph) {
    const mt_kahypar::HypernodeWeight perfectly_balanced_weight =
      std::ceil(partitioned_graph.totalWeight() / static_cast<double>(partitioned_graph.k()));
    double max_balance = partitioned_graph.partWeight(0) / static_cast<double>(perfectly_balanced_weight);
    for ( mt_kahypar::PartitionID i = 1; i < partitioned_graph.k(); ++i ) {
      max_balance = std::max(max_balance,
        partitioned_graph.partWeight(i) / static_cast<double>(perfectly_balanced_weight));
    }
    return max_balance - 1.0;
  }

  template<typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph partition(typename TypeTraits::Hypergraph& hypergraph,
                                                       Context& context) {
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    const bool is_graph = PartitionedHypergraph::TYPE == MULTILEVEL_GRAPH_PARTITIONING ||
                          PartitionedHypergraph::TYPE == N_LEVEL_HYPERGRAPH_PARTITIONING;
    if ( is_graph || context.partition.preset_type != PresetType::large_k ||
         PartitionedHypergraph::TYPE == LARGE_K_PARTITIONING ) {
      if ( lib::check_if_all_relavant_parameters_are_set(context) ) {
        mt_kahypar_hypergraph_t hg = utils::hypergraph_cast(hypergraph);
        if ( lib::check_compatibility(hg, lib::get_preset_c_type(context.partition.preset_type)) ) {
          context.partition.instance_type = lib::get_instance_type(hg);
          context.partition.partition_type = to_partition_c_type(
            context.partition.preset_type, context.partition.instance_type);
          lib::prepare_context(context);
          context.partition.num_vcycles = 0;
          return Partitioner<TypeTraits>::partition(hypergraph, context);
        } else {
          WARNING(lib::incompatibility_description(hg));
        }
      }
    } else {
      WARNING("You want to partition the hypergraph into a large number of blocks,"
        << "which is only possible when calling the function partitionIntoLargeK(...).");
    }
    return PartitionedHypergraph();
  }

  template<typename TypeTraits>
  void improve(typename TypeTraits::PartitionedHypergraph& partitioned_hg,
               Context& context,
               const size_t num_vcycles) {
    if ( lib::check_if_all_relavant_parameters_are_set(context) ) {
      mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
      if ( lib::check_compatibility(phg, lib::get_preset_c_type(context.partition.preset_type)) ) {
        context.partition.instance_type = lib::get_instance_type(phg);
        context.partition.partition_type = to_partition_c_type(
          context.partition.preset_type, context.partition.instance_type);
        lib::prepare_context(context);
        context.partition.num_vcycles = num_vcycles;
        Partitioner<TypeTraits>::partitionVCycle(partitioned_hg, context);
      } else {
        WARNING(lib::incompatibility_description(phg));
      }
    }
  }

}

PYBIND11_MODULE(mtkahypar, m) {

  // ####################### Enum Types #######################

  using mt_kahypar::FileFormat;
  py::enum_<FileFormat>(m, "FileFormat", py::module_local())
    .value("HMETIS", FileFormat::hMetis)
    .value("METIS", FileFormat::Metis);

  using mt_kahypar::PresetType;
  py::enum_<PresetType>(m, "PresetType", py::module_local())
    .value("DETERMINISTIC", PresetType::deterministic)
    .value("LARGE_K", PresetType::large_k)
    .value("DEFAULT", PresetType::default_preset)
    .value("DEFAULT_FLOWS", PresetType::default_flows);

  using mt_kahypar::Objective;
  py::enum_<Objective>(m, "Objective", py::module_local())
    .value("CUT", Objective::cut)
    .value("KM1", Objective::km1)
    .value("SOED", Objective::soed);

  // ####################### Initialize Thread Pool #######################

  m.def("initializeThreadPool", &initialize_thread_pool,
    "Initializes the thread pool with the given number of threads",
    py::arg("number of threads"));

  // ####################### Context #######################

  py::class_<Context>(m, "Context", py::module_local())
    .def(py::init<>())
    .def("loadPreset", [](Context& context, const PresetType preset) {
        switch ( preset ) {
          case PresetType::deterministic:
            context.load_deterministic_preset();
            break;
          case PresetType::large_k:
            context.load_large_k_preset();
            break;
          case PresetType::default_preset:
            context.load_default_preset();
            break;
          case PresetType::default_flows:
            context.load_default_flow_preset();
            break;
          default:
            LOG << "Preset type" << preset << "not supported!";
            break;
        }
      }, "Loads a preset for partitioning (DETERMINISTIC, SPEED or HIGH_QUALITY)",
      py::arg("preset type"))
    .def("loadConfigurationFile", [](Context& context, const std::string& config_file) {
        mt_kahypar::parseIniToContext(context, config_file);
      }, "Read partitioning configuration from file",
      py::arg("configuration file"))
    .def("setPartitioningParameters",
      [](Context& context,
         const PartitionID k,
         const double epsilon,
         const Objective objective,
         const size_t seed) {
           context.partition.k = k;
           context.partition.epsilon = epsilon;
           context.partition.objective = objective;
           context.partition.seed = seed;
         }, "Sets all required parameters for partitioning",
         py::arg("k"), py::arg("epsilon"), py::arg("objective function"), py::arg("seed"))
    .def_property("k",
      [](const Context& context) {
        return context.partition.k;
      }, [](Context& context, const PartitionID k) {
        context.partition.k = k;
      }, "Number of blocks in which the (hyper)graph should be partitioned into")
    .def_property("epsilon",
      [](const Context& context) {
        return context.partition.epsilon;
      }, [](Context& context, const double epsilon) {
        context.partition.epsilon = epsilon;
      }, "Allowed imbalance")
    .def_property("objective",
      [](const Context& context) {
        return context.partition.objective;
      }, [](Context& context, const Objective objective) {
        context.partition.objective = objective;
      }, "Sets the objective function for partitioning (CUT, KM1 or SOED)")
    .def_property("seed",
      [](const Context& context) {
        return context.partition.seed;
      }, [](Context& context, const int seed) {
        context.partition.seed = seed;
      }, "Seed for the random number generator")
    .def_property("num_vcycles",
      [](const Context& context) {
        return context.partition.num_vcycles;
      }, [](Context& context, const size_t num_vcycles) {
        context.partition.num_vcycles = num_vcycles;
      }, "Sets the number of V-cycles")
    .def_property("logging",
      [](const Context& context) {
        return context.partition.verbose_output;
      }, [](Context& context, const bool verbose_output) {
        context.partition.verbose_output = verbose_output;
      }, "Enable partitioning output")
    .def_property("max_block_weights",
      [](const Context& context) {
        return context.partition.max_part_weights;
      }, [](Context& context, std::vector<HypernodeWeight>& block_weights) {
        context.partition.use_individual_part_weights = true;
        context.partition.max_part_weights.assign(block_weights.size(), 0);
        for ( size_t block = 0; block < block_weights.size(); ++block ) {
          context.partition.max_part_weights[block] = block_weights[block];
        }
      }, "Maximum allowed weight for each block of the output partition")
    .def("outputConfiguration", [](const Context& context) {
        LOG << context;
      }, "Output partitioning configuration");

  // ####################### Graph #######################

  using Graph = ds::StaticGraph;
  using GraphFactory = typename Graph::Factory;
  py::class_<Graph>(m, "Graph")
    .def(py::init<>([](const HypernodeID num_nodes,
                       const HyperedgeID num_edges,
                       const vec<std::pair<HypernodeID,HypernodeID>>& edges) {
        return GraphFactory::construct_from_graph_edges(
          num_nodes, num_edges, edges, nullptr, nullptr, true);
      }), R"pbdoc(
Construct an unweighted graph.

:param num_nodes: Number of nodes
:param num_edges: Number of edges
:param edges: list of tuples containing all edges (e.g., [(0,1),(0,2),(1,3),...])
          )pbdoc",
      py::arg("num_nodes"),
      py::arg("num_edges"),
      py::arg("edges"))
    .def(py::init<>([](const HypernodeID num_nodes,
                       const HyperedgeID num_edges,
                       const vec<std::pair<HypernodeID,HypernodeID>>& edges,
                       const vec<HypernodeWeight>& node_weights,
                       const vec<HyperedgeWeight>& edge_weights) {
        return GraphFactory::construct_from_graph_edges(
          num_nodes, num_edges, edges, edge_weights.data(), node_weights.data(), true);
      }), R"pbdoc(
Construct a weighted graph.

:param num_nodes: Number of nodes
:param num_edges: Number of edges
:param edges: list of tuples containing all edges (e.g., [(0,1),(0,2),(1,3),...])
:param node_weights: Weights of all nodes
:param hyperedge_weights: Weights of all edges
          )pbdoc",
      py::arg("num_nodes"),
      py::arg("num_edges"),
      py::arg("edges"),
      py::arg("node_weights"),
      py::arg("edge_weights"))
    .def(py::init<>([](const std::string& file_name,
                      const FileFormat file_format) {
        return io::readInputFile<Graph>(file_name, file_format, true);
      }), "Reads a graph from a file (supported file formats are METIS and HMETIS)",
      py::arg("filename"), py::arg("format"))
    .def("numNodes", &Graph::initialNumNodes,
      "Number of nodes")
    .def("numEdges", [](Graph& graph) {
        return graph.initialNumEdges() / 2;
      }, "Number of undirected edges")
    .def("numDirectedEdges", &Graph::initialNumEdges,
      "Number of directed edges")
    .def("totalWeight", &Graph::totalWeight,
      "Total weight of all nodes")
    .def("nodeDegree", &Graph::nodeDegree,
      "Degree of node", py::arg("node"))
    .def("nodeWeight", &Graph::nodeWeight,
      "Weight of node", py::arg("node"))
    .def("edgeWeight", &Graph::edgeWeight,
      "Weight of edge", py::arg("edge"))
    .def("source", &Graph::edgeSource,
      "Source node of edge (e.g., (0,1) -> 0 is the source node)",
      py::arg("edge"))
    .def("target", &Graph::edgeTarget,
      "Target node of edge (e.g., (0,1) -> 1 is the target node)",
      py::arg("edge"))
    .def("doForAllNodes",
      [&](Graph& graph,
          const std::function<void(const HypernodeID&)>& f) {
        for ( const HypernodeID& hn : graph.nodes() ) {
          f(hn);
        }
      }, "Executes lambda expression for all nodes",
      py::arg("lambda"))
    .def("doForAllEdges",
      [&](Graph& graph,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : graph.edges() ) {
          f(he);
        }
      }, "Executes lambda expression for all edges",
      py::arg("lambda"))
    .def("doForAllIncidentEdges",
      [&](Graph& graph,
          const HypernodeID hn,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : graph.incidentEdges(hn) ) {
          f(he);
        }
      }, "Executes lambda expression for all incident edges of a node",
      py::arg("node"), py::arg("lambda"))
    .def("doForAllNeighbors",
      [&](Graph& graph,
          const HypernodeID hn,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : graph.incidentEdges(hn) ) {
          f(graph.edgeTarget(he));
        }
      }, "Executes lambda expression for all adjacent nodes of a node",
      py::arg("node"), py::arg("lambda"))
    .def("partition", &partition<StaticGraphTypeTraits>,
      "Partitions the graph with the parameters given in the corresponding context",
      py::arg("context"));

  // ####################### Hypergraph #######################

  using Hypergraph = ds::StaticHypergraph;
  using HypergraphFactory = typename Hypergraph::Factory;
  py::class_<Hypergraph>(m, "Hypergraph")
    .def(py::init<>([](const HypernodeID num_hypernodes,
                       const HyperedgeID num_hyperedges,
                       const vec<vec<HypernodeID>>& hyperedges) {
        return HypergraphFactory::construct(
          num_hypernodes, num_hyperedges, hyperedges, nullptr, nullptr, true);
      }), R"pbdoc(
Construct an unweighted hypergraph.

:param num_hypernodes: Number of nodes
:param num_hyperedges: Number of hyperedges
:param hyperedges: list containing all hyperedges (e.g., [[0,1],[0,2,3],...])
          )pbdoc",
      py::arg("num_hypernodes"),
      py::arg("num_hyperedges"),
      py::arg("hyperedges"))
    .def(py::init<>([](const HypernodeID num_hypernodes,
                       const HyperedgeID num_hyperedges,
                       const vec<vec<HypernodeID>>& hyperedges,
                       const vec<HypernodeWeight>& node_weights,
                       const vec<HyperedgeWeight>& hyperedge_weights) {
        return HypergraphFactory::construct(
          num_hypernodes, num_hyperedges, hyperedges,
          hyperedge_weights.data(), node_weights.data(), true);
      }), R"pbdoc(
Construct a weighted hypergraph.

:param num_hypernodes: Number of nodes
:param num_hyperedges: Number of hyperedges
:param hyperedges: List containing all hyperedges (e.g., [[0,1],[0,2,3],...])
:param node_weights: Weights of all hypernodes
:param hyperedge_weights: Weights of all hyperedges
          )pbdoc",
      py::arg("num_hypernodes"),
      py::arg("num_hyperedges"),
      py::arg("hyperedges"),
      py::arg("node_weights"),
      py::arg("hyperedge_weights"))
    .def(py::init<>([](const std::string& file_name,
                       const FileFormat file_format) {
        return io::readInputFile<Hypergraph>(file_name, file_format, true);
      }), "Reads a hypergraph from a file (supported file formats are METIS and HMETIS)",
      py::arg("filename"), py::arg("format"))
    .def("numNodes", &Hypergraph::initialNumNodes,
      "Number of nodes")
    .def("numEdges", &Hypergraph::initialNumEdges,
      "Number of hyperedges")
    .def("numPins", &Hypergraph::initialNumPins,
      "Number of pins")
    .def("totalWeight", &Hypergraph::totalWeight,
      "Total weight of all nodes")
    .def("nodeDegree", &Hypergraph::nodeDegree,
      "Degree of node", py::arg("node"))
    .def("nodeWeight", &Hypergraph::nodeWeight,
      "Weight of node", py::arg("node"))
    .def("edgeSize", &Hypergraph::edgeSize,
      "Size of hyperedge", py::arg("hyperedge"))
    .def("edgeWeight", &Hypergraph::edgeWeight,
      "Weight of hyperedge", py::arg("hyperedge"))
    .def("doForAllNodes",
      [&](Hypergraph& hypergraph,
          const std::function<void(const HypernodeID&)>& f) {
        for ( const HypernodeID& hn : hypergraph.nodes() ) {
          f(hn);
        }
      }, "Executes lambda expression for all nodes",
      py::arg("lambda"))
    .def("doForAllEdges",
      [&](Hypergraph& hypergraph,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : hypergraph.edges() ) {
          f(he);
        }
      }, "Executes lambda expression for all hyperedges",
      py::arg("lambda"))
    .def("doForAllIncidentEdges",
      [&](Hypergraph& hypergraph,
          const HypernodeID hn,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          f(he);
        }
      }, "Executes lambda expression for all incident hyperedges of a node",
      py::arg("node"), py::arg("lambda"))
    .def("doForAllPins",
      [&](Hypergraph& hypergraph,
          const HyperedgeID& he,
          const std::function<void(const HypernodeID&)>& f) {
        for ( const HyperedgeID& hn : hypergraph.pins(he) ) {
          f(hn);
        }
      }, "Executes lambda expression for all pins of a hyperedge",
      py::arg("hyperedge"), py::arg("lambda"))
    .def("partition", &partition<StaticHypergraphTypeTraits>,
      "Partitions the hypergraph with the parameters given in the corresponding context",
      py::arg("context"))
    .def("partitionIntoLargeK", &partition<LargeKHypergraphTypeTraits>,
      "Partitions the hypergraph into a large number of blocks with the parameters given in the corresponding context",
      py::arg("context"));

  // ####################### Partitioned Graph #######################

  using PartitionedGraph = typename StaticGraphTypeTraits::PartitionedHypergraph;
  py::class_<PartitionedGraph>(m, "PartitionedGraph")
    .def(py::init<>([](Graph& graph,
                       const PartitionID num_blocks,
                       const vec<PartitionID>& partition) {
        PartitionedGraph partitioned_graph(num_blocks, graph, parallel_tag_t { });
        partitioned_graph.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
          }
          partitioned_graph.setOnlyNodePart(hn, partition[hn]);
        });
        partitioned_graph.initializePartition();
        return partitioned_graph;
      }), R"pbdoc(
Construct a partitioned graph.

:param graph: graph object
:param num_blocks: number of block in which the graph should be partitioned into
:param partition: List of block IDs for each node
          )pbdoc",
      py::arg("graph"), py::arg("num_blocks"), py::arg("partition"))
    .def(py::init<>([](Graph& graph,
                       const PartitionID num_blocks,
                       const std::string& partition_file) {
        std::vector<PartitionID> partition;
        io::readPartitionFile(partition_file, partition);
        PartitionedGraph partitioned_graph(num_blocks, graph, mt_kahypar::parallel_tag_t { });
        partitioned_graph.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
          }
          partitioned_graph.setOnlyNodePart(hn, partition[hn]);
        });
        partitioned_graph.initializePartition();
        return partitioned_graph;
      }), R"pbdoc(
Construct a partitioned graph.

:param graph: graph object
:param num_blocks: number of block in which the graph should be partitioned into
:param partition_file: Partition file containing block IDs for each node
          )pbdoc",
      py::arg("graph"), py::arg("num_blocks"), py::arg("partition_file"))
    .def("numBlocks", &PartitionedGraph::k,
      "Number of blocks")
    .def("blockWeight", &PartitionedGraph::partWeight,
      "Weight of the corresponding block", py::arg("block"))
    .def("blockID", &PartitionedGraph::partID,
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("isIncidentToCutEdge", &PartitionedGraph::isBorderNode,
      "Returns true, if the corresponding node is incident to at least one cut edge",
      py::arg("node"))
    .def("numIncidentCutEdges", &PartitionedGraph::numIncidentCutHyperedges,
      "Number of incident cut edges of the corresponding node",
      py::arg("node"))
    .def("connectivity", &PartitionedGraph::connectivity,
      "Either one, if the corresponding edge is non-cut edge, or two if it is a cut edge",
      py::arg("edge"))
    .def("doForAllBlocksInEdge",
      [](PartitionedGraph& partitioned_graph,
         const HyperedgeID he,
         const std::function<void(const PartitionID&)>& f) {
        for ( const PartitionID& block : partitioned_graph.connectivitySet(he) ) {
          f(block);
        }
      }, "Executes lambda expression on blocks contained in the given edge",
      py::arg("edge"), py::arg("lambda"))
    .def("imbalance", [](PartitionedGraph& partitioned_graph) {
        return imbalance(partitioned_graph);
      }, "Computes the imbalance of the partition")
    .def("cut", [](PartitionedGraph& partitioned_graph) {
        return metrics::quality(partitioned_graph, Objective::cut);
      },
      "Computes the edge-cut metric of the partition")
    .def("writePartitionToFile", [](PartitionedGraph& partitioned_graph,
                                    const std::string& partition_file) {
        io::writePartitionFile(partitioned_graph, partition_file);
      }, "Writes the partition to a file",
      py::arg("partition_file"))
    .def("improve", &improve<StaticGraphTypeTraits>,
      "Improves the partition using the iterated multilevel cycle technique (V-cycles)",
      py::arg("context"), py::arg("num_vcycles"));

  // ####################### Partitioned Hypergraph #######################

  using PartitionedHypergraph = typename StaticHypergraphTypeTraits::PartitionedHypergraph;
  py::class_<PartitionedHypergraph>(m, "PartitionedHypergraph")
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const vec<PartitionID>& partition) {
        PartitionedHypergraph partitioned_hg(num_blocks, hypergraph, parallel_tag_t { });
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
          }
          partitioned_hg.setOnlyNodePart(hn, partition[hn]);
        });
        partitioned_hg.initializePartition();
        return partitioned_hg;
      }), R"pbdoc(
Construct a partitioned hypergraph.

:param hypergraph: hypergraph object
:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition: List of block IDs for each node
          )pbdoc",
      py::arg("hypergraph"), py::arg("num_blocks"), py::arg("partition"))
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const std::string& partition_file) {
        std::vector<PartitionID> partition;
        io::readPartitionFile(partition_file, partition);
        PartitionedHypergraph partitioned_hg(num_blocks, hypergraph, parallel_tag_t { });
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
          }
          partitioned_hg.setOnlyNodePart(hn, partition[hn]);
        });
        partitioned_hg.initializePartition();
        return partitioned_hg;
      }), R"pbdoc(
Construct a partitioned hypergraph.

:param hypergraph: hypergraph object
:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition_file: Partition file containing block IDs for each node
          )pbdoc",
      py::arg("hypergraph"), py::arg("num_blocks"), py::arg("partition_file"))
    .def("numBlocks", &PartitionedHypergraph::k,
      "Number of blocks")
    .def("blockWeight", &PartitionedHypergraph::partWeight,
      "Weight of the corresponding block", py::arg("block"))
    .def("blockID", &PartitionedHypergraph::partID,
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("isIncidentToCutEdge", &PartitionedHypergraph::isBorderNode,
      "Returns true, if the corresponding node is incident to at least one cut hyperedge",
      py::arg("node"))
    .def("numIncidentCutEdges", &PartitionedHypergraph::numIncidentCutHyperedges,
      "Number of incident cut hyperedges of the corresponding node",
      py::arg("node"))
    .def("numPinsInBlock", &PartitionedHypergraph::pinCountInPart,
      "Number of nodes part of the corresponding block in the given hyperedge",
      py::arg("hyperedge"), py::arg("block"))
    .def("connectivity", &PartitionedHypergraph::connectivity,
      "Number of distinct blocks to which the pins of corresponding hyperedge are assigned",
      py::arg("hyperedge"))
    .def("doForAllBlocksInEdge",
      [](PartitionedHypergraph& partitioned_hg,
         const HyperedgeID he,
         const std::function<void(const PartitionID&)>& f) {
        for ( const PartitionID& block : partitioned_hg.connectivitySet(he) ) {
          f(block);
        }
      }, "Executes lambda expression on blocks contained in the given hyperedge",
      py::arg("hyperedge"), py::arg("lambda"))
    .def("imbalance", [](PartitionedHypergraph& partitioned_hg) {
        return imbalance(partitioned_hg);
      }, "Computes the imbalance of the partition")
    .def("cut", [](PartitionedHypergraph& partitioned_hg) {
        return metrics::quality(partitioned_hg, Objective::cut);
      },
      "Computes the cut-net metric of the partition")
    .def("km1", [](PartitionedHypergraph& partitioned_hg) {
        return metrics::quality(partitioned_hg, Objective::km1);
      },
      "Computes the connectivity metric of the partition")
    .def("soed", [](PartitionedHypergraph& partitioned_hg) {
        return metrics::quality(partitioned_hg, Objective::soed);
      },
      "Computes the sum-of-external-degree metric of the partition")
    .def("writePartitionToFile", [](PartitionedHypergraph& partitioned_hg,
                                    const std::string& partition_file) {
        io::writePartitionFile(partitioned_hg, partition_file);
      }, "Writes the partition to a file",
      py::arg("partition_file"))
    .def("improve", &improve<StaticHypergraphTypeTraits>,
      "Improves the partition using the iterated multilevel cycle technique (V-cycles)",
      py::arg("context"), py::arg("num_vcycles"));

 // ####################### Partitioned Hypergraph #######################

  using SparsePartitionedHypergraph = typename LargeKHypergraphTypeTraits::PartitionedHypergraph;
  py::class_<SparsePartitionedHypergraph>(m, "SparsePartitionedHypergraph")
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const vec<PartitionID>& partition) {
        SparsePartitionedHypergraph partitioned_hg(num_blocks, hypergraph, parallel_tag_t { });
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
          }
          partitioned_hg.setOnlyNodePart(hn, partition[hn]);
        });
        partitioned_hg.initializePartition();
        return partitioned_hg;
      }), R"pbdoc(
Construct a partitioned hypergraph.

:param hypergraph: hypergraph object
:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition: List of block IDs for each node
          )pbdoc",
      py::arg("hypergraph"), py::arg("num_blocks"), py::arg("partition"))
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const std::string& partition_file) {
        std::vector<PartitionID> partition;
        io::readPartitionFile(partition_file, partition);
        SparsePartitionedHypergraph partitioned_hg(num_blocks, hypergraph, parallel_tag_t { });
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
          }
          partitioned_hg.setOnlyNodePart(hn, partition[hn]);
        });
        partitioned_hg.initializePartition();
        return partitioned_hg;
      }), R"pbdoc(
Construct a partitioned hypergraph.

:param hypergraph: hypergraph object
:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition_file: Partition file containing block IDs for each node
          )pbdoc",
      py::arg("hypergraph"), py::arg("num_blocks"), py::arg("partition_file"))
    .def("numBlocks", &SparsePartitionedHypergraph::k,
      "Number of blocks")
    .def("blockWeight", &SparsePartitionedHypergraph::partWeight,
      "Weight of the corresponding block", py::arg("block"))
    .def("blockID", &SparsePartitionedHypergraph::partID,
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("isIncidentToCutEdge", &SparsePartitionedHypergraph::isBorderNode,
      "Returns true, if the corresponding node is incident to at least one cut hyperedge",
      py::arg("node"))
    .def("numIncidentCutEdges", &SparsePartitionedHypergraph::numIncidentCutHyperedges,
      "Number of incident cut hyperedges of the corresponding node",
      py::arg("node"))
    .def("numPinsInBlock", &SparsePartitionedHypergraph::pinCountInPart,
      "Number of nodes part of the corresponding block in the given hyperedge",
      py::arg("hyperedge"), py::arg("block"))
    .def("connectivity", &SparsePartitionedHypergraph::connectivity,
      "Number of distinct blocks to which the pins of corresponding hyperedge are assigned",
      py::arg("hyperedge"))
    .def("doForAllBlocksInEdge",
      [](SparsePartitionedHypergraph& partitioned_hg,
         const HyperedgeID he,
         const std::function<void(const PartitionID&)>& f) {
        for ( const PartitionID& block : partitioned_hg.connectivitySet(he) ) {
          f(block);
        }
      }, "Executes lambda expression on blocks contained in the given hyperedge",
      py::arg("hyperedge"), py::arg("lambda"))
    .def("imbalance", [](SparsePartitionedHypergraph& partitioned_hg) {
        return imbalance(partitioned_hg);
      }, "Computes the imbalance of the partition")
    .def("cut", [](SparsePartitionedHypergraph& partitioned_hg) {
        return metrics::quality(partitioned_hg, Objective::cut);
      },
      "Computes the cut-net metric of the partition")
    .def("km1", [](SparsePartitionedHypergraph& partitioned_hg) {
        return metrics::quality(partitioned_hg, Objective::km1);
      },
      "Computes the connectivity metric of the partition")
    .def("soed", [](SparsePartitionedHypergraph& partitioned_hg) {
        return metrics::quality(partitioned_hg, Objective::soed);
      },
      "Computes the sum-of-external-degree metric of the partition")
    .def("writePartitionToFile", [](SparsePartitionedHypergraph& partitioned_hg,
                                    const std::string& partition_file) {
        io::writePartitionFile(partitioned_hg, partition_file);
      }, "Writes the partition to a file",
      py::arg("partition_file"))
    .def("improve", &improve<LargeKHypergraphTypeTraits>,
      "Improves the partition using the iterated multilevel cycle technique (V-cycles)",
      py::arg("context"), py::arg("num_vcycles"));


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

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
#include <sstream>

#ifndef KAHYPAR_DISABLE_HWLOC
  #include <hwloc.h>
#endif

#include "include/mtkahypartypes.h"
#include "include/helper_functions.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/partitioner.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/io/hypergraph_factory.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/io/presets.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/io/command_line_options.h"


namespace py = pybind11;
using namespace mt_kahypar;

namespace {
  void initialize(const size_t num_threads) {
    lib::initialize(num_threads, true);
  }

  void throw_if_not_compatible(mt_kahypar_partitioned_hypergraph_t partitioned_hg, PresetType preset) {
    if (!lib::check_compatibility(partitioned_hg, lib::get_preset_c_type(preset))) {
      throw UnsupportedOperationException(lib::incompatibility_description(partitioned_hg));
    }
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
  typename TypeTraits::PartitionedHypergraph createPartitionedHypergraph(typename TypeTraits::Hypergraph& hg,
                                                                         const PartitionID num_blocks,
                                                                         const std::vector<PartitionID>& partition) {
    typename TypeTraits::PartitionedHypergraph phg(num_blocks, hg, parallel_tag_t { });
    phg.doParallelForAllNodes([&](const HypernodeID& hn) {
      PartitionID block = partition[hn];
      if ( block < 0 || block >= num_blocks ) {
        std::stringstream ss;
        ss << "Invalid block ID for node " << hn << " ( block ID = " << partition[hn] << " )";
        throw InvalidInputException(ss.str());
      }
      phg.setOnlyNodePart(hn, block);
    });
    phg.initializePartition();
    return phg;
  }


  // ####################### Partitioning #######################

  template<typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph partitionImpl(typename TypeTraits::Hypergraph& hypergraph,
                                                           Context& context,
                                                           TargetGraph* target_graph) {
    using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
    const bool is_graph = PartitionedHypergraph::TYPE == MULTILEVEL_GRAPH_PARTITIONING ||
                          PartitionedHypergraph::TYPE == N_LEVEL_HYPERGRAPH_PARTITIONING;
    if ( !is_graph && context.partition.preset_type == PresetType::large_k &&
         PartitionedHypergraph::TYPE != LARGE_K_PARTITIONING ) {
      throw UnsupportedOperationException(
        "For hypergraphs, large k partitioning is only possible with the function partitionIntoLargeK(...).");
    }
    mt_kahypar_hypergraph_t hg = utils::hypergraph_cast(hypergraph);
    if ( !lib::check_compatibility(hg, lib::get_preset_c_type(context.partition.preset_type)) ) {
      throw UnsupportedOperationException(lib::incompatibility_description(hg));
    }
    if ( !lib::check_if_all_relavant_parameters_are_set(context) ) {
      throw InvalidInputException("Some required context parameter is not set.");
    }

    context.partition.instance_type = lib::get_instance_type(hg);
    context.partition.partition_type = to_partition_c_type(
      context.partition.preset_type, context.partition.instance_type);
    lib::prepare_context(context);
    context.partition.num_vcycles = 0;
    return Partitioner<TypeTraits>::partition(hypergraph, context, target_graph);
  }

  template<typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph partition(typename TypeTraits::Hypergraph& hypergraph,
                                                       const Context& context) {
    Context partition_context(context);
    return partitionImpl<TypeTraits>(hypergraph, partition_context, nullptr);
  }

  template<typename TypeTraits>
  typename TypeTraits::PartitionedHypergraph map(typename TypeTraits::Hypergraph& hypergraph,
                                                 const ds::StaticGraph& graph,
                                                 const Context& context) {
    TargetGraph target_graph(graph.copy(parallel_tag_t { }));
    Context partition_context(context);
    partition_context.partition.objective = Objective::steiner_tree;
    return partitionImpl<TypeTraits>(hypergraph, partition_context, &target_graph);
  }


  // ####################### V-Cycles #######################

  template<typename TypeTraits>
  void improveImpl(typename TypeTraits::PartitionedHypergraph& partitioned_hg,
                   Context& context,
                   const size_t num_vcycles,
                   TargetGraph* target_graph) {
    mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
    if ( !lib::check_compatibility(phg, lib::get_preset_c_type(context.partition.preset_type)) ) {
      throw UnsupportedOperationException(lib::incompatibility_description(phg));
    }
    if ( !lib::check_if_all_relavant_parameters_are_set(context) ) {
      throw InvalidInputException("Some required context parameter is not set.");
    }

    context.partition.instance_type = lib::get_instance_type(phg);
    context.partition.partition_type = to_partition_c_type(
      context.partition.preset_type, context.partition.instance_type);
    lib::prepare_context(context);
    context.partition.num_vcycles = num_vcycles;
    Partitioner<TypeTraits>::partitionVCycle(partitioned_hg, context, target_graph);
  }

  template<typename TypeTraits>
  void improve(typename TypeTraits::PartitionedHypergraph& partitioned_hg,
               const Context& context,
               const size_t num_vcycles) {
    Context partition_context(context);
    improveImpl<TypeTraits>(partitioned_hg, partition_context, num_vcycles, nullptr);
  }

  template<typename TypeTraits>
  void improveMapping(typename TypeTraits::PartitionedHypergraph& partitioned_hg,
                      const ds::StaticGraph& graph,
                      const Context& context,
                      const size_t num_vcycles) {
    TargetGraph target_graph(graph.copy(parallel_tag_t { }));
    Context partition_context(context);
    partition_context.partition.objective = Objective::steiner_tree;
    improveImpl<TypeTraits>(partitioned_hg, partition_context, num_vcycles, &target_graph);
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
    .value("QUALITY", PresetType::quality);

  using mt_kahypar::Objective;
  py::enum_<Objective>(m, "Objective", py::module_local())
    .value("CUT", Objective::cut)
    .value("KM1", Objective::km1)
    .value("SOED", Objective::soed);

  // ####################### Exceptions #######################

  py::register_exception<InvalidInputException>(m, "InvalidInputError", PyExc_ValueError);
  py::register_exception<InvalidParameterException>(m, "InvalidParameterError", PyExc_ValueError);
  py::register_exception<UnsupportedOperationException>(m, "UnsupportedOperationError");
  py::register_exception<SystemException>(m, "SystemError");

  // ####################### Initialize Thread Pool #######################

  m.def("initialize", [&](const size_t num_threads) {
      initialize(num_threads);
    }, "General initialization. Initializes the thread pool with the given number of threads",
    py::arg("number of threads"));

  // ####################### Initialize Random Number Generator #######################

  m.def("set_seed", [&](const int seed) {
      mt_kahypar::utils::Randomize::instance().setSeed(seed);
    }, "Initializes the random number generator with the given seed",
    py::arg("seed"));

  // ####################### Context #######################

  py::class_<Context>(m, "Context", py::module_local())
    .def(py::init<>([](const PresetType preset) {
        Context context;
        auto preset_option_list = loadPreset(preset);
        mt_kahypar::presetToContext(context, preset_option_list, true);
        return context;
      }))
    .def(py::init<>([](const std::string& config_file) {
        Context context;
        mt_kahypar::parseIniToContext(context, config_file, true);
        return context;
      }))
    .def("set_partitioning_parameters",
      [](Context& context,
         const PartitionID k,
         const double epsilon,
         const Objective objective) {
           context.partition.k = k;
           context.partition.epsilon = epsilon;
           context.partition.objective = objective;
         }, "Sets all required parameters for partitioning",
         py::arg("k"), py::arg("epsilon"), py::arg("objective function"))
    .def("set_mapping_parameters",
      [](Context& context,
         const PartitionID k,
         const double epsilon) {
           context.partition.k = k;
           context.partition.epsilon = epsilon;
         }, "Sets all required parameters for mapping to a target graph",
         py::arg("k"), py::arg("epsilon"))
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
    .def("print_configuration", [](const Context& context) {
        LOG << context;
      }, "Print partitioning configuration");

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
    .def("num_nodes", &Graph::initialNumNodes,
      "Number of nodes")
    .def("num_edges", [](Graph& graph) {
        return graph.initialNumEdges() / 2;
      }, "Number of undirected edges")
    .def("num_directed_edges", &Graph::initialNumEdges,
      "Number of directed edges")
    .def("total_weight", &Graph::totalWeight,
      "Total weight of all nodes")
    .def("node_degree", &Graph::nodeDegree,
      "Degree of node", py::arg("node"))
    .def("node_weight", &Graph::nodeWeight,
      "Weight of node", py::arg("node"))
    .def("is_fixed", &Graph::isFixed,
      "Returns whether or not the corresponding node is a fixed vertex",
      py::arg("node"))
    .def("fixed_vertex_block", [&](const Graph& graph,
                                 const HypernodeID hn) {
        return graph.isFixed(hn) ? graph.fixedVertexBlock(hn) : kInvalidPartition;
      }, "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("add_fixed_vertices", [&](Graph& graph,
                                 const vec<PartitionID>& fixed_vertices,
                                 const PartitionID num_blocks) {
        mt_kahypar_hypergraph_t gr = utils::hypergraph_cast(graph);
        io::addFixedVertices(gr, fixed_vertices.data(), num_blocks);
      }, R"pbdoc(
Adds the fixed vertices specified in the array to the graph. The array must contain
n entries (n = number of nodes). Each entry contains either the fixed vertex block of the
corresponding node or -1 if the node is not fixed.
          )pbdoc", py::arg("fixed_vertices"), py::arg("num_blocks"))
    .def("add_fixed_vertices_from_file", [&](Graph& graph,
                                         const std::string& fixed_vertex_file,
                                         const PartitionID num_blocks) {
        mt_kahypar_hypergraph_t gr = utils::hypergraph_cast(graph);
        io::addFixedVerticesFromFile(gr, fixed_vertex_file, num_blocks);
      }, R"pbdoc(
Adds the fixed vertices specified in the fixed vertex file to the graph. The file must contain
n lines (n = number of nodes). Each line contains either the fixed vertex block of the
corresponding node or -1 if the node is not fixed.
          )pbdoc", py::arg("fixed_vertex_file"), py::arg("num_blocks"))
    .def("remove_fixed_vertices", [&](Graph& graph) {
        mt_kahypar_hypergraph_t gr = utils::hypergraph_cast(graph);
        io::removeFixedVertices(gr);
    }, "Removes all fixed vertices from the hypergraph")
    .def("edge_weight", &Graph::edgeWeight,
      "Weight of edge", py::arg("edge"))
    .def("source", &Graph::edgeSource,
      "Source node of edge (e.g., (0,1) -> 0 is the source node)",
      py::arg("edge"))
    .def("target", &Graph::edgeTarget,
      "Target node of edge (e.g., (0,1) -> 1 is the target node)",
      py::arg("edge"))
    .def("do_for_all_nodes",
      [&](Graph& graph,
          const std::function<void(const HypernodeID&)>& f) {
        for ( const HypernodeID& hn : graph.nodes() ) {
          f(hn);
        }
      }, "Executes lambda expression for all nodes",
      py::arg("lambda"))
    .def("do_for_all_edges",
      [&](Graph& graph,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : graph.edges() ) {
          f(he);
        }
      }, "Executes lambda expression for all edges",
      py::arg("lambda"))
    .def("do_for_all_incident_edges",
      [&](Graph& graph,
          const HypernodeID hn,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : graph.incidentEdges(hn) ) {
          f(he);
        }
      }, "Executes lambda expression for all incident edges of a node",
      py::arg("node"), py::arg("lambda"))
    .def("do_for_all_neighbors",
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
      py::arg("context"))
    .def("map_onto_graph", &map<StaticGraphTypeTraits>,
      R"pbdoc(
  Maps a (hyper)graph onto a target graph with the configuration specified in the partitioning context.
  The number of blocks of the output mapping/partition is the same as the number of nodes in the target graph
  (each node of the target graph represents a block). The objective is to minimize the total weight of
  all Steiner trees spanned by the (hyper)edges on the target graph. A Steiner tree is a tree with minimal weight
  that spans a subset of the nodes (in our case the hyperedges) on the target graph. This objective function
  is able to acurately model wire-lengths in VLSI design or communication costs in a distributed system where some
  processors do not communicate directly with each other or different speeds.
          )pbdoc", py::arg("target_graph"), py::arg("context"));

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
    .def("num_nodes", &Hypergraph::initialNumNodes,
      "Number of nodes")
    .def("num_edges", &Hypergraph::initialNumEdges,
      "Number of hyperedges")
    .def("num_pins", &Hypergraph::initialNumPins,
      "Number of pins")
    .def("total_weight", &Hypergraph::totalWeight,
      "Total weight of all nodes")
    .def("node_degree", &Hypergraph::nodeDegree,
      "Degree of node", py::arg("node"))
    .def("node_weight", &Hypergraph::nodeWeight,
      "Weight of node", py::arg("node"))
    .def("is_fixed", &Hypergraph::isFixed,
      "Returns whether or not the corresponding node is a fixed vertex",
      py::arg("node"))
    .def("fixed_vertex_block", [&](const Hypergraph& hypergraph,
                                 const HypernodeID hn) {
        return hypergraph.isFixed(hn) ? hypergraph.fixedVertexBlock(hn) : kInvalidPartition;
      }, "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("add_fixed_vertices", [&](Hypergraph& hypergraph,
                                 const vec<PartitionID>& fixed_vertices,
                                 const PartitionID num_blocks) {
        mt_kahypar_hypergraph_t hg = utils::hypergraph_cast(hypergraph);
        io::addFixedVertices(hg, fixed_vertices.data(), num_blocks);
      }, R"pbdoc(
Adds the fixed vertices specified in the array to the hypergraph. The array must contain
n entries (n = number of nodes). Each entry contains either the fixed vertex block of the
corresponding node or -1 if the node is not fixed.
          )pbdoc", py::arg("fixed_vertices"), py::arg("num_blocks"))
    .def("add_fixed_vertices_from_file", [&](Hypergraph& hypergraph,
                                         const std::string& fixed_vertex_file,
                                         const PartitionID num_blocks) {
        mt_kahypar_hypergraph_t hg = utils::hypergraph_cast(hypergraph);
        io::addFixedVerticesFromFile(hg, fixed_vertex_file, num_blocks);
      }, R"pbdoc(
Adds the fixed vertices specified in the fixed vertex file to the hypergraph. The file must contain
n lines (n = number of nodes). Each line contains either the fixed vertex block of the
corresponding node or -1 if the node is not fixed.
          )pbdoc", py::arg("fixed_vertex_file"), py::arg("num_blocks"))
    .def("remove_fixed_vertices", [&](Hypergraph& hypergraph) {
        mt_kahypar_hypergraph_t hg = utils::hypergraph_cast(hypergraph);
        io::removeFixedVertices(hg);
    }, "Removes all fixed vertices from the hypergraph")
    .def("edge_size", &Hypergraph::edgeSize,
      "Size of hyperedge", py::arg("hyperedge"))
    .def("edge_weight", &Hypergraph::edgeWeight,
      "Weight of hyperedge", py::arg("hyperedge"))
    .def("do_for_all_nodes",
      [&](Hypergraph& hypergraph,
          const std::function<void(const HypernodeID&)>& f) {
        for ( const HypernodeID& hn : hypergraph.nodes() ) {
          f(hn);
        }
      }, "Executes lambda expression for all nodes",
      py::arg("lambda"))
    .def("do_for_all_edges",
      [&](Hypergraph& hypergraph,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : hypergraph.edges() ) {
          f(he);
        }
      }, "Executes lambda expression for all hyperedges",
      py::arg("lambda"))
    .def("do_for_all_incident_edges",
      [&](Hypergraph& hypergraph,
          const HypernodeID hn,
          const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          f(he);
        }
      }, "Executes lambda expression for all incident hyperedges of a node",
      py::arg("node"), py::arg("lambda"))
    .def("do_for_all_pins",
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
    .def("partition_into_large_k", &partition<LargeKHypergraphTypeTraits>,
      "Partitions the hypergraph into a large number of blocks with the parameters given in the corresponding context",
      py::arg("context"))
    .def("map_onto_graph", &map<StaticHypergraphTypeTraits>,
      R"pbdoc(
  Maps a (hyper)graph onto a target graph with the configuration specified in the partitioning context.
  The number of blocks of the output mapping/partition is the same as the number of nodes in the target graph
  (each node of the target graph represents a block). The objective is to minimize the total weight of
  all Steiner trees spanned by the (hyper)edges on the target graph. A Steiner tree is a tree with minimal weight
  that spans a subset of the nodes (in our case the hyperedges) on the target graph. This objective function
  is able to acurately model wire-lengths in VLSI design or communication costs in a distributed system where some
  processors do not communicate directly with each other or different speeds.
          )pbdoc", py::arg("target_graph"), py::arg("context"));

  // ####################### Partitioned Graph #######################

  using PartitionedGraph = typename StaticGraphTypeTraits::PartitionedHypergraph;
  py::class_<PartitionedGraph>(m, "PartitionedGraph")
    .def(py::init<>([](Graph& graph,
                       const PartitionID num_blocks,
                       const std::vector<PartitionID>& partition) {
        return createPartitionedHypergraph<StaticGraphTypeTraits>(graph, num_blocks, partition);
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
        return createPartitionedHypergraph<StaticGraphTypeTraits>(graph, num_blocks, partition);
      }), R"pbdoc(
Construct a partitioned graph.

:param graph: graph object
:param num_blocks: number of block in which the graph should be partitioned into
:param partition_file: Partition file containing block IDs for each node
          )pbdoc",
      py::arg("graph"), py::arg("num_blocks"), py::arg("partition_file"))
    .def("num_blocks", &PartitionedGraph::k,
      "Number of blocks")
    .def("block_weight", &PartitionedGraph::partWeight,
      "Weight of the corresponding block", py::arg("block"))
    .def("block_id", &PartitionedGraph::partID,
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("is_fixed", &PartitionedGraph::isFixed,
      "Returns whether or not the corresponding node is a fixed vertex",
      py::arg("node"))
    .def("fixed_vertex_block", [&](const PartitionedGraph& graph,
                                 const HypernodeID hn) {
        return graph.isFixed(hn) ? graph.fixedVertexBlock(hn) : kInvalidPartition;
      }, "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("is_incident_to_cut_edge", &PartitionedGraph::isBorderNode,
      "Returns true, if the corresponding node is incident to at least one cut edge",
      py::arg("node"))
    .def("num_incident_cut_edges", &PartitionedGraph::numIncidentCutHyperedges,
      "Number of incident cut edges of the corresponding node",
      py::arg("node"))
    .def("connectivity", &PartitionedGraph::connectivity,
      "Either one, if the corresponding edge is non-cut edge, or two if it is a cut edge",
      py::arg("edge"))
    .def("do_for_all_blocks_in_edge",
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
    .def("steiner_tree", [](PartitionedGraph& partitioned_graph,
                            Graph& graph) {
        TargetGraph target_graph(graph.copy(parallel_tag_t { }));
        target_graph.precomputeDistances(4);
        partitioned_graph.setTargetGraph(&target_graph);
        return metrics::quality(partitioned_graph, Objective::steiner_tree);
      }, "Computes the steiner tree metric of the mapping",
      py::arg("target_graph"))
    .def("write_partition_to_file", [](PartitionedGraph& partitioned_graph,
                                    const std::string& partition_file) {
        io::writePartitionFile(partitioned_graph, partition_file);
      }, "Writes the partition to a file",
      py::arg("partition_file"))
    .def("improve_partition", &improve<StaticGraphTypeTraits>,
      "Improves the partition using the iterated multilevel cycle technique (V-cycles)",
      py::arg("context"), py::arg("num_vcycles"))
    .def("improve_mapping", &improveMapping<StaticGraphTypeTraits>,
      "Improves a mapping onto a graph using the iterated multilevel cycle technique (V-cycles)",
      py::arg("target_graph"), py::arg("context"), py::arg("num_vcycles"));

  // ####################### Partitioned Hypergraph #######################

  using PartitionedHypergraph = typename StaticHypergraphTypeTraits::PartitionedHypergraph;
  py::class_<PartitionedHypergraph>(m, "PartitionedHypergraph")
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const std::vector<PartitionID>& partition) {
        return createPartitionedHypergraph<StaticHypergraphTypeTraits>(hypergraph, num_blocks, partition);
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
        return createPartitionedHypergraph<StaticHypergraphTypeTraits>(hypergraph, num_blocks, partition);
      }), R"pbdoc(
Construct a partitioned hypergraph.

:param hypergraph: hypergraph object
:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition_file: Partition file containing block IDs for each node
          )pbdoc",
      py::arg("hypergraph"), py::arg("num_blocks"), py::arg("partition_file"))
    .def("num_blocks", &PartitionedHypergraph::k,
      "Number of blocks")
    .def("block_weight", &PartitionedHypergraph::partWeight,
      "Weight of the corresponding block", py::arg("block"))
    .def("block_id", &PartitionedHypergraph::partID,
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("is_fixed", &PartitionedHypergraph::isFixed,
      "Returns whether or not the corresponding node is a fixed vertex",
      py::arg("node"))
    .def("fixed_vertex_block", [&](const PartitionedHypergraph& hypergraph,
                                 const HypernodeID hn) {
        return hypergraph.isFixed(hn) ? hypergraph.fixedVertexBlock(hn) : kInvalidPartition;
      }, "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("is_incident_to_cut_edge", &PartitionedHypergraph::isBorderNode,
      "Returns true, if the corresponding node is incident to at least one cut hyperedge",
      py::arg("node"))
    .def("num_incident_cut_edges", &PartitionedHypergraph::numIncidentCutHyperedges,
      "Number of incident cut hyperedges of the corresponding node",
      py::arg("node"))
    .def("num_pins_in_block", &PartitionedHypergraph::pinCountInPart,
      "Number of nodes part of the corresponding block in the given hyperedge",
      py::arg("hyperedge"), py::arg("block"))
    .def("connectivity", &PartitionedHypergraph::connectivity,
      "Number of distinct blocks to which the pins of corresponding hyperedge are assigned",
      py::arg("hyperedge"))
    .def("do_for_all_blocks_in_edge",
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
    .def("steiner_tree", [](PartitionedHypergraph& partitioned_hg,
                            Graph& graph) {
        TargetGraph target_graph(graph.copy(parallel_tag_t { }));
        target_graph.precomputeDistances(4);
        partitioned_hg.setTargetGraph(&target_graph);
        return metrics::quality(partitioned_hg, Objective::steiner_tree);
      }, "Computes the steiner tree metric of the mapping",
      py::arg("target_graph"))
    .def("write_partition_to_file", [](PartitionedHypergraph& partitioned_hg,
                                    const std::string& partition_file) {
        io::writePartitionFile(partitioned_hg, partition_file);
      }, "Writes the partition to a file",
      py::arg("partition_file"))
    .def("improve_partition", &improve<StaticHypergraphTypeTraits>,
      "Improves the partition using the iterated multilevel cycle technique (V-cycles)",
      py::arg("context"), py::arg("num_vcycles"))
    .def("improve_mapping", &improveMapping<StaticHypergraphTypeTraits>,
      "Improves a mapping onto a graph using the iterated multilevel cycle technique (V-cycles)",
      py::arg("target_graph"), py::arg("context"), py::arg("num_vcycles"));

 // ####################### Partitioned Hypergraph #######################

  using SparsePartitionedHypergraph = typename LargeKHypergraphTypeTraits::PartitionedHypergraph;
  py::class_<SparsePartitionedHypergraph>(m, "SparsePartitionedHypergraph")
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const std::vector<PartitionID>& partition) {
        return createPartitionedHypergraph<LargeKHypergraphTypeTraits>(hypergraph, num_blocks, partition);
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
        return createPartitionedHypergraph<LargeKHypergraphTypeTraits>(hypergraph, num_blocks, partition);
      }), R"pbdoc(
Construct a partitioned hypergraph.

:param hypergraph: hypergraph object
:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition_file: Partition file containing block IDs for each node
          )pbdoc",
      py::arg("hypergraph"), py::arg("num_blocks"), py::arg("partition_file"))
    .def("num_blocks", &SparsePartitionedHypergraph::k,
      "Number of blocks")
    .def("block_weight", &SparsePartitionedHypergraph::partWeight,
      "Weight of the corresponding block", py::arg("block"))
    .def("block_id", &SparsePartitionedHypergraph::partID,
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("is_fixed", &SparsePartitionedHypergraph::isFixed,
      "Returns whether or not the corresponding node is a fixed vertex",
      py::arg("node"))
    .def("fixed_vertex_block", [&](const SparsePartitionedHypergraph& hypergraph,
                                 const HypernodeID hn) {
        return hypergraph.isFixed(hn) ? hypergraph.fixedVertexBlock(hn) : kInvalidPartition;
      }, "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("is_incident_to_cut_edge", &SparsePartitionedHypergraph::isBorderNode,
      "Returns true, if the corresponding node is incident to at least one cut hyperedge",
      py::arg("node"))
    .def("num_incident_cut_edges", &SparsePartitionedHypergraph::numIncidentCutHyperedges,
      "Number of incident cut hyperedges of the corresponding node",
      py::arg("node"))
    .def("num_pins_in_block", &SparsePartitionedHypergraph::pinCountInPart,
      "Number of nodes part of the corresponding block in the given hyperedge",
      py::arg("hyperedge"), py::arg("block"))
    .def("connectivity", &SparsePartitionedHypergraph::connectivity,
      "Number of distinct blocks to which the pins of corresponding hyperedge are assigned",
      py::arg("hyperedge"))
    .def("do_for_all_blocks_in_edge",
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
    .def("write_partition_to_file", [](SparsePartitionedHypergraph& partitioned_hg,
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

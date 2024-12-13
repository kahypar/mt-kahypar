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

#include <atomic>
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
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/delete.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/utils/randomize.h"


namespace py = pybind11;
using namespace mt_kahypar;

namespace {
  // Token type for initialization
  struct Initializer {};

  // Graph types
  struct mt_kahypar_py_graph_t: public mt_kahypar_hypergraph_t {};

  struct mt_kahypar_py_target_graph_t: public mt_kahypar_py_graph_t {};

  void initialize(const size_t num_threads, const bool print_warnings) {
    static std::atomic_bool is_initialized = false;
    bool expected = false;
    if (is_initialized.compare_exchange_strong(expected, true)) {
      lib::initialize(num_threads, true, print_warnings);
    } else if (print_warnings) {
      WARNING("Mt-KaHyPar is already initialized");
    }
  }

  const ds::StaticGraph& target_graph_cast(mt_kahypar_py_target_graph_t target_graph) {
    ASSERT(target_graph.type == STATIC_GRAPH);
    return utils::cast_const<ds::StaticGraph>(target_graph);
  }

  // Deleters
  struct HypergraphDeleter {
    void operator()(mt_kahypar_hypergraph_t* hg) {
      utils::delete_hypergraph(*hg);
      std::default_delete<mt_kahypar_hypergraph_t>{}(hg);
    }
  };

  struct PartitionedHypergraphDeleter {
    void operator()(mt_kahypar_partitioned_hypergraph_t* phg) {
      utils::delete_partitioned_hypergraph(*phg);
      std::default_delete<mt_kahypar_partitioned_hypergraph_t>{}(phg);
    }
  };
}


PYBIND11_MODULE(mtkahypar, m) {
  using StaticGraphFactory = typename ds::StaticGraph::Factory;

  // forward declare all classes so that the name lookup for docu and error messages works correctly
  auto initializer_class = py::class_<Initializer>(m, "Initializer");

  auto context_class = py::class_<Context>(m, "Context");

  auto hg_class = py::class_<mt_kahypar_hypergraph_t,
    std::unique_ptr<mt_kahypar_hypergraph_t, HypergraphDeleter>>(m, "Hypergraph");

  auto graph_class = py::class_<mt_kahypar_py_graph_t, mt_kahypar_hypergraph_t,
    std::unique_ptr<mt_kahypar_py_graph_t, HypergraphDeleter>>(m, "Graph");

  py::class_<mt_kahypar_py_target_graph_t, mt_kahypar_py_graph_t,
    std::unique_ptr<mt_kahypar_py_target_graph_t, HypergraphDeleter>>(m, "TargetGraph");
  
  auto phg_class = py::class_<mt_kahypar_partitioned_hypergraph_t,
    std::unique_ptr<mt_kahypar_partitioned_hypergraph_t, PartitionedHypergraphDeleter>>(m, "PartitionedHypergraph");


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
    .value("QUALITY", PresetType::quality)
    .value("HIGHEST_QUALITY", PresetType::highest_quality);

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

  // ####################### Initialize Thread Pool and RNG #######################

  m.def("initialize", [&](const size_t num_threads, const bool print_warnings) {
      initialize(num_threads, print_warnings);
      return Initializer{};
    }, "Initializes Mt-KaHyPar with the given number of threads.",
    py::arg("number of threads"), py::arg("print_warnings") = true);

  m.def("set_seed", [&](const int seed) {
      mt_kahypar::utils::Randomize::instance().setSeed(seed);
    }, "Initializes the random number generator with the given seed",
    py::arg("seed"));


  // ####################### The Initializer #######################

  initializer_class
    .def("context_from_preset",
      [](Initializer&, const PresetType preset) {
        return lib::context_from_preset(preset);
      }, "Creates a context from the given preset.",
      py::arg("preset"))
    .def("context_from_file",
      [](Initializer&, const std::string& config_file) {
        return lib::context_from_file(config_file.c_str());
      }, "Creates a context from a configuration file.",
      py::arg("config_file"))
    .def("create_hypergraph",
      [](Initializer&,
         const Context& context,
         const HypernodeID num_hypernodes,
         const HyperedgeID num_hyperedges,
         const vec<vec<HypernodeID>>& hyperedges) {
        return lib::create_hypergraph(
          context, num_hypernodes, num_hyperedges, hyperedges, nullptr, nullptr);
      }, R"pbdoc(
Construct an unweighted hypergraph.

:param context: the partitioning context
:param num_hypernodes: Number of nodes
:param num_hyperedges: Number of hyperedges
:param hyperedges: list containing all hyperedges (e.g., [[0,1],[0,2,3],...])
          )pbdoc",
      py::arg("context"),
      py::arg("num_hypernodes"),
      py::arg("num_hyperedges"),
      py::arg("hyperedges"))
    // Note: Using default arguments, i.e., nullable raw pointers, for the weights does not seem to work here.
    // The overload resolution fails due to a bug in pybind (seems related to https://github.com/pybind/pybind11/issues/4956)
    .def("create_hypergraph",
      [](Initializer&,
         const Context& context,
         const HypernodeID num_hypernodes,
         const HyperedgeID num_hyperedges,
         const vec<vec<HypernodeID>>& hyperedges,
         const vec<HypernodeWeight>& node_weights,
         const vec<HyperedgeWeight>& hyperedge_weights) {
        return lib::create_hypergraph(
          context, num_hypernodes, num_hyperedges, hyperedges, hyperedge_weights.data(), node_weights.data());
      }, R"pbdoc(
Construct a weighted hypergraph.

:param context: the partitioning context
:param num_hypernodes: Number of nodes
:param num_hyperedges: Number of hyperedges
:param hyperedges: List containing all hyperedges (e.g., [[0,1],[0,2,3],...])
:param node_weights: Weights of all hypernodes
:param hyperedge_weights: Weights of all hyperedges
          )pbdoc",
      py::arg("context"),
      py::arg("num_hypernodes"),
      py::arg("num_hyperedges"),
      py::arg("hyperedges"),
      py::arg("node_weights"),
      py::arg("hyperedge_weights"))
    .def("hypergraph_from_file",
      [](Initializer&,
         const std::string& file_name,
         const Context& context,
         const FileFormat file_format) {
        const bool stable_construction = context.partition.preset_type == PresetType::deterministic ? true : false;
        return io::readInputFile(file_name, context.partition.preset_type, InstanceType::hypergraph, file_format, stable_construction);
      }, "Reads a hypergraph from a file (supported file formats are METIS and HMETIS)",
      py::arg("filename"), py::arg("context"), py::arg("format"))
    .def("create_graph",
      [](Initializer&,
         const Context& context,
         const HypernodeID num_nodes,
         const HyperedgeID num_edges,
         const vec<std::pair<HypernodeID,HypernodeID>>& edges) {
        return mt_kahypar_py_graph_t{lib::create_graph(
          context, num_nodes, num_edges, edges, nullptr, nullptr)};
      }, R"pbdoc(
Construct an unweighted graph.

:param context: the partitioning context
:param num_nodes: Number of nodes
:param num_edges: Number of edges
:param edges: list of tuples containing all edges (e.g., [(0,1),(0,2),(1,3),...])
          )pbdoc",
      py::arg("context"),
      py::arg("num_nodes"),
      py::arg("num_edges"),
      py::arg("edges"))
    .def("create_graph",
      [](Initializer&,
         const Context& context,
         const HypernodeID num_nodes,
         const HyperedgeID num_edges,
         const vec<std::pair<HypernodeID,HypernodeID>>& edges,
         const vec<HypernodeWeight>& node_weights,
         const vec<HyperedgeWeight>& edge_weights) {
        return mt_kahypar_py_graph_t{lib::create_graph(
          context, num_nodes, num_edges, edges, edge_weights.data(), node_weights.data())};
      }, R"pbdoc(
Construct a weighted graph.

:param context: the partitioning context
:param num_nodes: Number of nodes
:param num_edges: Number of edges
:param edges: list of tuples containing all edges (e.g., [(0,1),(0,2),(1,3),...])
:param node_weights: Weights of all nodes
:param edge_weights: Weights of all edges
          )pbdoc",
      py::arg("context"),
      py::arg("num_nodes"),
      py::arg("num_edges"),
      py::arg("edges"),
      py::arg("node_weights"),
      py::arg("edge_weights"))
    .def("graph_from_file",
      [](Initializer&,
         const std::string& file_name,
         const Context& context,
         const FileFormat file_format) {
        const bool stable_construction = context.partition.preset_type == PresetType::deterministic ? true : false;
        return mt_kahypar_py_graph_t{io::readInputFile(
          file_name, context.partition.preset_type, InstanceType::graph, file_format, stable_construction)};
      }, "Reads a graph from a file (supported file formats are METIS and HMETIS)",
      py::arg("filename"), py::arg("context"), py::arg("format"))
    .def("create_target_graph",
      [](Initializer&,
         const Context& context,
         const HypernodeID num_nodes,
         const HyperedgeID num_edges,
         const vec<std::pair<HypernodeID,HypernodeID>>& edges,
         const vec<HyperedgeWeight>& edge_weights) {
        // TODO: check input lengths
        unused(context);
        return mt_kahypar_py_target_graph_t{
          reinterpret_cast<mt_kahypar_hypergraph_s*>(new ds::StaticGraph(
            StaticGraphFactory::construct_from_graph_edges(num_nodes, num_edges,
              edges, edge_weights.data(), nullptr, true))), STATIC_GRAPH };
      }, R"pbdoc(
Construct a target graph.

:param context: the partitioning context
:param num_nodes: Number of nodes
:param num_edges: Number of edges
:param edges: list of tuples containing all edges (e.g., [(0,1),(0,2),(1,3),...])
:param edge_weights: Weights of all edges
          )pbdoc",
      py::arg("context"),
      py::arg("num_nodes"),
      py::arg("num_edges"),
      py::arg("edges"),
      py::arg("edge_weights"))
    .def("target_graph_from_file",
      [](Initializer&,
         const std::string& file_name,
         const Context& context,
         const FileFormat file_format) {
        unused(context);
        return mt_kahypar_py_target_graph_t{
          reinterpret_cast<mt_kahypar_hypergraph_s*>(new ds::StaticGraph(
            io::readInputFile<ds::StaticGraph>(file_name, file_format, true))),
            STATIC_GRAPH };
      }, "Reads a target graph from a file (supported file formats are METIS and HMETIS)",
      py::arg("filename"), py::arg("context"), py::arg("format"));


  // ####################### Context #######################

  context_class
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
        if (k != context.partition.k) {
          context.partition.k = k;
          context.partition.use_individual_part_weights = false;
        }
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
    .def("set_individual_target_block_weights",
      [](Context& context, std::vector<HypernodeWeight>& block_weights) {
        if (static_cast<PartitionID>(block_weights.size()) != context.partition.k) {
          throw InvalidParameterException("Number of block weights must be equal to k");
        }
        lib::set_individual_block_weights(context, block_weights.size(), block_weights.data());
      }, "Set individual maximum allowed weight for each block of the output partition",
      py::arg("block_weights"))
    .def("compute_max_block_weights",
      [](Context& context, HypernodeWeight total_weight) {
        context.setupPartWeights(total_weight);
        return context.partition.max_part_weights;
      }, R"pbdoc(
Compute maximum allowed block weights for an input with the given total weight (depends on epsilon).
If indiviual target weights are set, these are returned instead.

:param total_weight: Total weight of input hypergraph
          )pbdoc",
      py::arg("total_weight"))
    .def("compute_perfect_balance_block_weights",
      [](Context& context, HypernodeWeight total_weight) {
        context.setupPartWeights(total_weight);
        return context.partition.perfect_balance_part_weights;
      }, R"pbdoc(
Compute block weights that represent perfect balance for an input with the given total weight.
If indiviual target weights are set, these are used for the calculation instead.

:param total_weight: Total weight of input hypergraph
          )pbdoc",
      py::arg("total_weight"))
    .def("print_configuration", [](const Context& context) {
        LOG << context;
      }, "Print partitioning configuration");

  // ####################### Hypergraph and Graph #######################

  // TODO: preset compatibility check
  hg_class
    .def("num_nodes", &lib::num_nodes<true>, "Number of nodes")
    .def("num_edges", &lib::num_edges<true>, "Number of hyperedges")
    .def("num_pins", &lib::num_pins<true>, "Number of pins")
    .def("total_weight", &lib::total_weight<true>, "Total weight of all nodes")
    .def("node_degree", &lib::node_degree<true>, "Degree of node", py::arg("node"))
    .def("node_weight", &lib::node_weight<true>, "Weight of node", py::arg("node"))
    .def("is_fixed", &lib::is_fixed<true>,
      "Returns whether or not the corresponding node is a fixed vertex", py::arg("node"))
    .def("fixed_vertex_block", &lib::fixed_vertex_block<true>,
      "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("add_fixed_vertices", [&](mt_kahypar_hypergraph_t hypergraph,
                                   const vec<PartitionID>& fixed_vertices,
                                   const PartitionID num_blocks) {
        io::addFixedVertices(hypergraph, fixed_vertices.data(), num_blocks);
      }, R"pbdoc(
Adds the fixed vertices specified in the array to the (hyper)graph. The array must contain
n entries (n = number of nodes). Each entry contains either the fixed vertex block of the
corresponding node or -1 if the node is not fixed.
          )pbdoc", py::arg("fixed_vertices"), py::arg("num_blocks"))
    .def("add_fixed_vertices_from_file", [&](mt_kahypar_hypergraph_t hypergraph,
                                             const std::string& fixed_vertex_file,
                                             const PartitionID num_blocks) {
        io::addFixedVerticesFromFile(hypergraph, fixed_vertex_file, num_blocks);
      }, R"pbdoc(
Adds the fixed vertices specified in the fixed vertex file to the (hyper)graph. The file must contain
n lines (n = number of nodes). Each line contains either the fixed vertex block of the
corresponding node or -1 if the node is not fixed.
          )pbdoc", py::arg("fixed_vertex_file"), py::arg("num_blocks"))
    .def("remove_fixed_vertices", [&](mt_kahypar_hypergraph_t hypergraph) {
        io::removeFixedVertices(hypergraph);
    }, "Removes all fixed vertices from the hypergraph")
    .def("edge_size", &lib::edge_size<true>,
         "Size of hyperedge", py::arg("hyperedge"))
    .def("edge_weight", &lib::edge_weight<true>,
         "Weight of hyperedge", py::arg("hyperedge"))
    // .def("do_for_all_nodes",
    //   [&](HypergraphT& hypergraph,
    //       const std::function<void(const HypernodeID&)>& f) {
    //     for ( const HypernodeID& hn : hypergraph.nodes() ) {
    //       f(hn);
    //     }
    //   }, "Executes lambda expression for all nodes",
    //   py::arg("lambda"))
    // .def("do_for_all_edges",
    //   [&](HypergraphT& hypergraph,
    //       const std::function<void(const HyperedgeID&)>& f) {
    //     for ( const HyperedgeID& he : hypergraph.edges() ) {
    //       f(he);
    //     }
    //   }, "Executes lambda expression for all hyperedges",
    //   py::arg("lambda"))
    // .def("do_for_all_incident_edges",
    //   [&](HypergraphT& hypergraph,
    //       const HypernodeID hn,
    //       const std::function<void(const HyperedgeID&)>& f) {
    //     for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
    //       f(he);
    //     }
    //   }, "Executes lambda expression for all incident hyperedges of a node",
    //   py::arg("node"), py::arg("lambda"))
    // .def("do_for_all_pins",
    //   [&](HypergraphT& hypergraph,
    //       const HyperedgeID& he,
    //       const std::function<void(const HypernodeID&)>& f) {
    //     for ( const HyperedgeID& hn : hypergraph.pins(he) ) {
    //       f(hn);
    //     }
    //   }, "Executes lambda expression for all pins of a hyperedge",
    //   py::arg("hyperedge"), py::arg("lambda"))
    .def("partition",
      [&](mt_kahypar_hypergraph_t hypergraph, const Context& context) {
        return lib::partition(hypergraph, context);
      }, "Partitions the hypergraph with the parameters given in the partitioning context",
      py::arg("context"))
    .def("map_onto_graph",
      [&](mt_kahypar_hypergraph_t hypergraph, mt_kahypar_py_target_graph_t target_graph, const Context& context) {
        return lib::map(hypergraph, target_graph_cast(target_graph), context);
      },
      R"pbdoc(
  Maps a (hyper)graph onto a target graph with the parameters given in the partitioning context.
  The number of blocks of the output mapping/partition is the same as the number of nodes in the target graph
  (each node of the target graph represents a block). The objective is to minimize the total weight of
  all Steiner trees spanned by the (hyper)edges on the target graph. A Steiner tree is a tree with minimal weight
  that spans a subset of the nodes (in our case the hyperedges) on the target graph. This objective function
  is able to acurately model wire-lengths in VLSI design or communication costs in a distributed system where some
  processors do not communicate directly with each other or different speeds.
          )pbdoc", py::arg("target_graph"), py::arg("context"))
  .def("create_partitioned_hypergraph",
    [&](mt_kahypar_hypergraph_t hypergraph,
        const Context& context,
        const PartitionID num_blocks,
        const std::vector<PartitionID>& partition) {
      auto result = lib::create_partitioned_hypergraph(hypergraph, context, num_blocks, partition.data());
      if (result.partitioned_hg == nullptr) {
        throw UnsupportedOperationException("Input is not a valid hypergraph!");
      }
      return result;
    }, R"pbdoc(
Construct a partitioned hypergraph from this hypergraph.

:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition: list of block IDs for each node
        )pbdoc",
    py::arg("num_blocks"), py::arg("context"),py::arg("partition"),
    // prevent hypergraph from being freed while the PHG is still alive
    py::keep_alive<0, 1>())
  .def("partitioned_hypergraph_from_file",
    [&](mt_kahypar_hypergraph_t hypergraph,
        const Context& context,
        const PartitionID num_blocks,
        const std::string& partition_file) {
      std::vector<PartitionID> partition;
      io::readPartitionFile(partition_file, lib::num_nodes<true>(hypergraph), partition);
      auto result = lib::create_partitioned_hypergraph(hypergraph, context, num_blocks, partition.data());
      if (result.partitioned_hg == nullptr) {
        throw UnsupportedOperationException("Input is not a valid hypergraph!");
      }
      return result;
    }, R"pbdoc(
Construct a partitioned hypergraph from this hypergraph.

:param num_blocks: number of block in which the hypergraph should be partitioned into
:param partition_file: partition file containing block IDs for each node
        )pbdoc",
    py::arg("num_blocks"), py::arg("context"), py::arg("partition_file"),
    // prevent hypergraph from being freed while the PHG is still alive
    py::keep_alive<0, 1>());


  graph_class
    .def("num_directed_edges",
      [&](mt_kahypar_py_graph_t g) {
        return lib::switch_graph<HyperedgeID, true>(g, [](const auto& graph) {
          return graph.initialNumEdges();
        });
      }, "Number of directed edges (equal to num_edges)")
    .def("num_undirected_edges",
      [&](mt_kahypar_py_graph_t g) {
        return lib::switch_graph<HyperedgeID, true>(g, [](const auto& graph) {
          return graph.initialNumEdges() / 2;
        });
      }, "Number of undirected edges (equal to num_edges / 2)")
    .def("edge_source",
      [&](mt_kahypar_py_graph_t g, HyperedgeID edge) {
        return lib::switch_graph<HypernodeID, true>(g, [=](const auto& graph) {
          return graph.edgeSource(edge);
        });
      }, "Source node of edge (e.g., (0,1) -> 0 is the source node)", py::arg("edge"))
    .def("edge_target",
      [&](mt_kahypar_py_graph_t g, HyperedgeID edge) {
        return lib::switch_graph<HypernodeID, true>(g, [=](const auto& graph) {
          return graph.edgeTarget(edge);
        });
      }, "Target node of edge (e.g., (0,1) -> 1 is the target node)", py::arg("edge"));
    // .def("do_for_all_neighbors",
    //   [&](Graph& graph,
    //       const HypernodeID hn,
    //       const std::function<void(const HyperedgeID&)>& f) {
    //     for ( const HyperedgeID& he : graph.incidentEdges(hn) ) {
    //       f(graph.edgeTarget(he));
    //     }
    //   }, "Executes lambda expression for all adjacent nodes of a node",
    //   py::arg("node"), py::arg("lambda"));

  // ####################### Partitioned Hypergraph #######################

  phg_class
    .def("num_blocks",
      [&](mt_kahypar_partitioned_hypergraph_t p) {
        return lib::switch_phg_throwing<PartitionID>(p, [](const auto& phg) {
          return phg.k();
        });
      },
      "Number of blocks")
    .def("block_weight",
      [&](mt_kahypar_partitioned_hypergraph_t p, const PartitionID block) {
        return lib::switch_phg_throwing<HypernodeWeight>(p, [=](const auto& phg) {
          return phg.partWeight(block);
        });
      },
      "Weight of the corresponding block", py::arg("block"))
    .def("block_id",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HypernodeID node) {
        return lib::switch_phg_throwing<HypernodeWeight>(p, [=](const auto& phg) {
          return phg.partID(node);
        });
      },
      "Block to which the corresponding node is assigned", py::arg("node"))
    .def("is_fixed",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HypernodeID node) {
        return lib::switch_phg_throwing<bool>(p, [=](const auto& phg) {
          return phg.isFixed(node);
        });
      },
      "Returns whether or not the corresponding node is a fixed vertex", py::arg("node"))
    .def("fixed_vertex_block",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HypernodeID node) {
        return lib::switch_phg_throwing<PartitionID>(p, [=](const auto& phg) {
          return phg.isFixed(node) ? phg.fixedVertexBlock(node) : -1;
        });
      },
      "Block to which the node is fixed (-1 if not fixed)", py::arg("node"))
    .def("is_incident_to_cut_edge",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HypernodeID node) {
        return lib::switch_phg_throwing<bool>(p, [=](const auto& phg) {
          return phg.isBorderNode(node);
        });
      },
      "Returns true, if the corresponding node is incident to at least one cut hyperedge", py::arg("node"))
    .def("num_incident_cut_edges",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HypernodeID node) {
        return lib::switch_phg_throwing<HyperedgeID>(p, [=](const auto& phg) {
          return phg.numIncidentCutHyperedges(node);
        });
      },
      "Number of incident cut hyperedges of the corresponding node", py::arg("node"))
    .def("connectivity",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HyperedgeID he) {
        return lib::switch_phg_throwing<PartitionID>(p, [=](const auto& phg) {
          return phg.connectivity(he);
        });
      },
      "Number of distinct blocks to which the pins of corresponding hyperedge are assigned", py::arg("hyperedge"))
    .def("num_pins_in_block",
      [&](mt_kahypar_partitioned_hypergraph_t p, const HyperedgeID he, const PartitionID block_id) {
        return lib::switch_phg_throwing<HypernodeID>(p, [=](const auto& phg) {
          return phg.pinCountInPart(he, block_id);
        });
      },
      "Number of nodes part of the corresponding block in the given hyperedge", py::arg("hyperedge"), py::arg("block_id"))
    .def("imbalance",
      [&](mt_kahypar_partitioned_hypergraph_t p, const Context& context) {
        return lib::switch_phg_throwing<PartitionID>(p, [&](const auto& phg) {
          return lib::imbalance(phg, context);
        });
      },
      "Computes the imbalance of the partition", py::arg("context"))
    .def("cut",
      [&](mt_kahypar_partitioned_hypergraph_t p) {
        return lib::switch_phg_throwing<PartitionID>(p, [&](const auto& phg) {
          return metrics::quality(phg, Objective::cut);
        });
      },
      "Computes the cut-net metric of the partition")
    .def("km1",
      [&](mt_kahypar_partitioned_hypergraph_t p) {
        return lib::switch_phg_throwing<PartitionID>(p, [&](const auto& phg) {
          return metrics::quality(phg, Objective::km1);
        });
      },
      "Computes the connectivity metric of the partition")
    .def("soed",
      [&](mt_kahypar_partitioned_hypergraph_t p) {
        return lib::switch_phg_throwing<PartitionID>(p, [&](const auto& phg) {
          return metrics::quality(phg, Objective::soed);
        });
      },
      "Computes the sum-of-external-degree metric of the partition")
    .def("steiner_tree",
      [&](mt_kahypar_partitioned_hypergraph_t p, mt_kahypar_py_target_graph_t graph) {
        return lib::switch_phg_throwing<PartitionID>(p, [&](auto& phg) {
          TargetGraph target_graph(target_graph_cast(graph).copy());
          target_graph.precomputeDistances(4);
          phg.setTargetGraph(&target_graph);
          return metrics::quality(phg, Objective::steiner_tree);
        });
      },
      "Computes the sum-of-external-degree metric of the partition", py::arg("target_graph"))
    .def("write_partition_to_file",
      [&](mt_kahypar_partitioned_hypergraph_t p, const std::string& partition_file) {
        lib::switch_phg_throwing(p, [&](auto& phg) {
          io::writePartitionFile(phg, partition_file);
        });
      },
      "Writes the partition to a file", py::arg("partition_file"))
    .def("improve_partition", &lib::improve,
      "Improves the partition using the iterated multilevel cycle technique (V-cycles)",
      py::arg("context"), py::arg("num_vcycles"))
    .def("improve_mapping",
      [&](mt_kahypar_partitioned_hypergraph_t phg, mt_kahypar_py_target_graph_t target_graph, const Context& context, size_t num_vcycles) {
        lib::improveMapping(phg, target_graph_cast(target_graph), context, num_vcycles);
      },
      "Improves a mapping onto a graph using the iterated multilevel cycle technique (V-cycles)",
      py::arg("target_graph"), py::arg("context"), py::arg("num_vcycles"));
    // .def("do_for_all_blocks_in_edge",
    //   [](PartitionedHypergraphT& partitioned_hg,
    //     const HyperedgeID he,
    //     const std::function<void(const PartitionID&)>& f) {
    //     for ( const PartitionID& block : partitioned_hg.connectivitySet(he) ) {
    //       f(block);
    //     }
    //   }, "Executes lambda expression on blocks contained in the given hyperedge",
    //   py::arg("hyperedge"), py::arg("lambda"))


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

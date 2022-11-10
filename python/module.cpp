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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_io.h"

namespace py = pybind11;

namespace {
  template<typename T>
  using vec = mt_kahypar::parallel::scalable_vector<T>;

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
}

PYBIND11_MODULE(mtkahypar, m) {

  // ####################### Enum Types #######################

  using mt_kahypar::FileFormat;
  py::enum_<FileFormat>(m, "FileFormat")
    .value("HMETIS", FileFormat::hMetis)
    .value("METIS", FileFormat::Metis);

  using mt_kahypar::PresetType;
  py::enum_<PresetType>(m, "PresetType")
    .value("DETERMINISTIC", PresetType::deterministic)
    .value("SPEED", PresetType::default_preset)
    .value("HIGH_QUALITY", PresetType::default_flows);

  using mt_kahypar::Objective;
  py::enum_<Objective>(m, "Objective")
    .value("CUT", Objective::cut)
    .value("KM1", Objective::km1);

  // ####################### Initialize Thread Pool #######################

  m.def("initializeThreadPool", &initialize_thread_pool,
    "Initializes the thread pool with the given number of threads",
    py::arg("number of threads"));

  // ####################### Hypergraph #######################

  using mt_kahypar::HypernodeID;
  using mt_kahypar::HyperedgeID;
  using mt_kahypar::PartitionID;
  using mt_kahypar::HypernodeWeight;
  using mt_kahypar::HyperedgeWeight;
  using mt_kahypar::Hypergraph;

  py::class_<Hypergraph>(m, "Hypergraph")
    .def(py::init<>())
    .def("construct", [](Hypergraph& hypergraph,
                         const HypernodeID num_hypernodes,
                         const HyperedgeID num_hyperedges,
                         const vec<vec<HypernodeID>>& hyperedges) {
        hypergraph = mt_kahypar::HypergraphFactory::construct(
          num_hypernodes, num_hyperedges, hyperedges);
      }, R"pbdoc(
Construct an unweighted hypergraph.

:param num_hypernodes: Number of nodes
:param num_hyperedges: Number of hyperedges
:param hyperedges: list containing all hyperedges (e.g., [[0,1],[0,2,3],...])
          )pbdoc",
      py::arg("num_hypernodes"),
      py::arg("num_hyperedges"),
      py::arg("hyperedges"))
    .def("construct", [](Hypergraph& hypergraph,
                         const HypernodeID num_hypernodes,
                         const HyperedgeID num_hyperedges,
                         const vec<vec<HypernodeID>>& hyperedges,
                         const vec<HypernodeWeight>& node_weights,
                         const vec<HyperedgeWeight>& hyperedge_weights) {
        hypergraph = mt_kahypar::HypergraphFactory::construct(
          num_hypernodes, num_hyperedges, hyperedges,
          hyperedge_weights.data(), node_weights.data());
      }, R"pbdoc(
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
    .def("construct", [](Hypergraph& hypergraph,
                         const std::string& file_name,
                         const FileFormat file_format) {
        hypergraph = mt_kahypar::io::readInputFile(file_name, file_format, true);
      }, "Reads a hypergraph from a file (supported file formats are METIS and HMETIS)",
      py::arg("path to hypergraph file"), py::arg("file format"))
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
    .def("doForAllNodes", [&](Hypergraph& hypergraph,
                              const std::function<void(const HypernodeID&)>& f) {
        for ( const HypernodeID& hn : hypergraph.nodes() ) {
          f(hn);
        }
      }, "Executes lambda expression for all nodes",
      py::arg("lambda expression"))
    .def("doForAllEdges", [&](Hypergraph& hypergraph,
                              const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : hypergraph.edges() ) {
          f(he);
        }
      }, "Executes lambda expression for all edges",
      py::arg("lambda expression"))
    .def("doForAllIncidentEdges", [&](Hypergraph& hypergraph,
                                      const HypernodeID hn,
                                      const std::function<void(const HyperedgeID&)>& f) {
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          f(he);
        }
      }, "Executes lambda expression for all incident edges of a node",
      py::arg("node"), py::arg("lambda expression"))
    .def("doForAllPins", [&](Hypergraph& hypergraph,
                             const HyperedgeID& he,
                             const std::function<void(const HypernodeID&)>& f) {
        for ( const HyperedgeID& hn : hypergraph.pins(he) ) {
          f(hn);
        }
      }, "Executes lambda expression for all pins of a hyperedge",
      py::arg("hyperedge"), py::arg("lambda expression"));

  // ####################### Setup Context #######################

  using mt_kahypar::Context;
  py::class_<Context>(m, "Context")
    .def(py::init<>())
    .def("loadPreset", [](Context& context, const PresetType preset) {
        switch ( preset ) {
          case PresetType::deterministic:
            context.load_deterministic_preset();
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
      [](Context& context, const PartitionID k, const double epsilon,
         const Objective objective, const size_t seed) {
           context.partition.k = k;
           context.partition.epsilon = epsilon;
           context.partition.objective = objective;
           context.partition.seed = seed;
         }, "Sets all required parameters for partitioning",
         py::arg("k"), py::arg("epsilon"), py::arg("objective function"), py::arg("seed"))
    .def("setK", [](Context& context, const PartitionID k) {
        context.partition.k = k;
      }, "Number of blocks in which the (hyper)graph should be partitioned into",
      py::arg("k"))
    .def("setEpsilon", [](Context& context, const double epsilon) {
        context.partition.epsilon = epsilon;
      }, "Allowed imbalance",
      py::arg("epsilon"))
    .def("setIndividualBlockWeights",
      [](Context& context, const std::vector<HypernodeWeight>& block_weights) {
        context.partition.use_individual_part_weights = true;
        context.partition.max_part_weights.assign(block_weights.size(), 0);
        for ( size_t block = 0; block < block_weights.size(); ++block ) {
          context.partition.max_part_weights[block] = block_weights[block];
        }
      }, "Assigns each block of the partition an individual maximum allowed block weight",
      py::arg("individual target block weights"))
    .def("setObjective", [](Context& context, const Objective objective) {
        context.partition.objective = objective;
      }, "Sets the objective function for partitioning (either CUT or KM1)",
      py::arg("objective function"))
    .def("setSeed", [](Context& context, const int seed) {
        context.partition.seed = seed;
      }, "Seed for the random number generator",
      py::arg("seed"))
    .def("setNumberOfVCycles", [](Context& context, const size_t num_vcycles) {
        context.partition.num_vcycles = num_vcycles;
      }, "Sets the number of V-cycles",
      py::arg("number of vcycles"))
    .def("enableLogging", [](Context& context, const bool verbose) {
        context.partition.verbose_output = verbose;
      }, "Enable partitioning output",
      py::arg("bool"))
    .def("outputConfiguration", [](const Context& context) {
        LOG << context;
      }, "Output partitioning configuration");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

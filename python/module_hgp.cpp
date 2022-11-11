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
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/partitioner.h"
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

  double imbalance(const mt_kahypar::PartitionedHypergraph& partitioned_hg) {
    const mt_kahypar::HypernodeWeight perfectly_balanced_weight =
      std::ceil(partitioned_hg.totalWeight() / static_cast<double>(partitioned_hg.k()));
    double max_balance = partitioned_hg.partWeight(0) / static_cast<double>(perfectly_balanced_weight);
    for ( mt_kahypar::PartitionID i = 1; i < partitioned_hg.k(); ++i ) {
      max_balance = std::max(max_balance,
        partitioned_hg.partWeight(i) / static_cast<double>(perfectly_balanced_weight));
    }
    return max_balance - 1.0;
  }

  void prepare_context(mt_kahypar::Context& context) {
    context.partition.mode = mt_kahypar::Mode::direct;
    context.shared_memory.num_threads = mt_kahypar::TBBInitializer::instance().total_number_of_threads();
    context.utility_id = mt_kahypar::utils::Utilities::instance().registerNewUtilityObjects();
    mt_kahypar::utils::Randomize::instance().setSeed(context.partition.seed);

    context.partition.perfect_balance_part_weights.clear();
    if ( !context.partition.use_individual_part_weights ) {
      context.partition.max_part_weights.clear();
    }

    if ( context.partition.objective == mt_kahypar::Objective::cut &&
         context.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1 ) {
      context.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut;
    }

    if ( context.partition.objective == mt_kahypar::Objective::cut &&
         context.initial_partitioning.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1 ) {
      context.initial_partitioning.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut;
    }

    if ( context.partition.objective == mt_kahypar::Objective::km1 &&
         context.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut ) {
      context.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1;
    }

    if ( context.partition.objective == mt_kahypar::Objective::km1 &&
         context.initial_partitioning.refinement.label_propagation.algorithm ==
         mt_kahypar::LabelPropagationAlgorithm::label_propagation_cut ) {
      context.initial_partitioning.refinement.label_propagation.algorithm =
        mt_kahypar::LabelPropagationAlgorithm::label_propagation_km1;
    }
  }

  mt_kahypar::PartitionedHypergraph partition(mt_kahypar::Hypergraph& hypergraph, mt_kahypar::Context& context) {
    prepare_context(context);
    return mt_kahypar::partition(hypergraph, context);
  }

  void improve_partition(mt_kahypar::PartitionedHypergraph& partitioned_hg,
                         mt_kahypar::Context& context,
                         const size_t num_vcycles) {
    prepare_context(context);
    context.partition.num_vcycles = num_vcycles;
    mt_kahypar::partitionVCycle(partitioned_hg, context);
  }
}

PYBIND11_MODULE(mtkahyparhgp, m) {

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
  using mt_kahypar::HypernodeWeight;
  using mt_kahypar::HyperedgeWeight;
  using mt_kahypar::Hypergraph;
  py::class_<Hypergraph>(m, "Hypergraph")
    .def(py::init<>([](const HypernodeID num_hypernodes,
                       const HyperedgeID num_hyperedges,
                       const vec<vec<HypernodeID>>& hyperedges) {
        return mt_kahypar::HypergraphFactory::construct(
          num_hypernodes, num_hyperedges, hyperedges);
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
        return mt_kahypar::HypergraphFactory::construct(
          num_hypernodes, num_hyperedges, hyperedges,
          hyperedge_weights.data(), node_weights.data());
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
        return mt_kahypar::io::readInputFile(file_name, file_format, true);
      }), "Reads a hypergraph from a file (supported file formats are METIS and HMETIS)",
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

  // ####################### Partitioned Hypergraph #######################

  using mt_kahypar::PartitionID;
  using mt_kahypar::PartitionedHypergraph;
  py::class_<PartitionedHypergraph>(m, "PartitionedHypergraph")
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const vec<PartitionID>& partition) {
        PartitionedHypergraph partitioned_hg(num_blocks, hypergraph, mt_kahypar::parallel_tag_t { });
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERROR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
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
      py::arg("hypergraph"),
      py::arg("number of blocks"),
      py::arg("partition"))
    .def(py::init<>([](Hypergraph& hypergraph,
                       const PartitionID num_blocks,
                       const std::string& partition_file) {
        std::vector<PartitionID> partition;
        mt_kahypar::io::readPartitionFile(partition_file, partition);
        PartitionedHypergraph partitioned_hg(num_blocks, hypergraph, mt_kahypar::parallel_tag_t { });
        partitioned_hg.doParallelForAllNodes([&](const HypernodeID& hn) {
          if ( partition[hn] < 0 || partition[hn] >= num_blocks ) {
            ERROR("Invalid block ID for node" << hn << "( block ID =" << partition[hn] << ")");
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
      py::arg("hypergraph"),
      py::arg("number of blocks"),
      py::arg("partition file"))
    .def("numBlocks", &PartitionedHypergraph::k,
      "Number of blocks")
    .def("blockWeight", &PartitionedHypergraph::partWeight,
      "Weight of all nodes in corresponding block", py::arg("block"))
    .def("blockID", &PartitionedHypergraph::partID,
      "Block ID of node", py::arg("node"))
    .def("isIncidentToCutEdge", &PartitionedHypergraph::isBorderNode,
      "Returns true, if the corresponding node is incident to a cut hyperedge",
      py::arg("node"))
    .def("numIncidentCutEdges", &PartitionedHypergraph::numIncidentCutHyperedges,
      "Number of incident cut hyperedges of the given node",
      py::arg("node"))
    .def("numPinsInBlock", &PartitionedHypergraph::pinCountInPart,
      "Number of pins contained in the given block of a hyperedge",
      py::arg("hyperedge"), py::arg("block"))
    .def("connectivity", &PartitionedHypergraph::connectivity,
      "Number of blocks contained in the given hyperedge",
      py::arg("hyperedge"))
    .def("doForAllBlocksInEdge", [](PartitionedHypergraph& partitioned_hg,
                                    const HyperedgeID he,
                                    const std::function<void(const PartitionID&)>& f) {
        for ( const PartitionID& block : partitioned_hg.connectivitySet(he) ) {
          f(block);
        }
      }, "Executes lambda expression on blocks contained in the given hyperedge",
      py::arg("hyperedge"), py::arg("lambda expression"))
    .def("imbalance", [](PartitionedHypergraph& partitioned_hg) {
        return imbalance(partitioned_hg);
      }, "Computes the imbalance of the partition")
    .def("cut", [](PartitionedHypergraph& partitioned_hg) {
        return mt_kahypar::metrics::hyperedgeCut(partitioned_hg);
      },
      "Computes the cut-net metric for the partitioned hypergraph")
    .def("km1", [](PartitionedHypergraph& partitioned_hg) {
        return mt_kahypar::metrics::km1(partitioned_hg);
      },
      "Computes the connectivity metric for the partitioned hypergraph")
    .def("soed", [](PartitionedHypergraph& partitioned_hg) {
        return mt_kahypar::metrics::soed(partitioned_hg);
      },
      "Computes the sum-of-external-degree metric for the partitioned hypergraph")
    .def("writePartitionToFile", [](PartitionedHypergraph& partitioned_hg,
                                    const std::string& partition_file) {
        mt_kahypar::io::writePartitionFile(partitioned_hg, partition_file);
      }, "Writes the partition to a file",
      py::arg("target partition file"));

  // ####################### Partitioning #######################

  m.def(
    "partition", &partition,
    "Compute a k-way partition of the hypergraph",
    py::arg("hypergraph"), py::arg("context"));

  m.def(
    "improvePartition", &improve_partition,
    "Improves a k-way partition (using the V-cycle technique)",
    py::arg("partitioned hypergraph"), py::arg("context"), py::arg("number of V-cycles"));

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

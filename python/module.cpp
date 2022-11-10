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

#include <string>
#include <vector>
#include <iostream>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/hypergraph_io.h"

namespace py = pybind11;


PYBIND11_MODULE(mtkahypar, m) {
  using mt_kahypar::HypernodeID;
  using mt_kahypar::HyperedgeID;
  using mt_kahypar::PartitionID;
  using mt_kahypar::HypernodeWeight;
  using mt_kahypar::HyperedgeWeight;
  using mt_kahypar::HypernodeWeight;



  // ####################### Setup Context #######################

  using mt_kahypar::PresetType;
  py::enum_<PresetType>(m, "PresetType")
    .value("DETERMINISTIC", PresetType::deterministic)
    .value("SPEED", PresetType::default_preset)
    .value("HIGH_QUALITY", PresetType::default_flows);

  using mt_kahypar::Objective;
  py::enum_<Objective>(m, "Objective")
    .value("CUT", Objective::cut)
    .value("KM1", Objective::km1);

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

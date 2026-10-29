<h1 align="center">DynPart - Dynamic Graph and Hypergraph Partitioner</h1>

License|Zenodo
:--:|:--:
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)|[![DOI](https://zenodo.org/badge/205879380.svg)](https://zenodo.org/badge/latestdoi/205879380)

Table of Contents
-----------

   * [DynPart](#dynpart)
   * [Getting Started](#getting-started)
   * [Running DynPart](#running-dynpart)
   * [Input Format](#input-format)
   * [Reproducing Results](#reproducing-results)
   * [Licensing](#licensing)


DynPart
-----------
DynPart is a dynamic hypergraph partitioning algorithm that can efficiently maintain a high-quality balanced partition of a hypergraph under dynamic changes. Changes include the insertion and deletion of hyperedges, vertices and pins. DynPart computes a initial partition of the hypergraph using the Mt-KaHyPar framework. It utilizes incremental refinement, localized FM-based refinement and V-Cycles to move vertices between blocks on the fly and maintain a high-quality partition.

### Mt-KaHyPar
DynPart is based on the [Mt-KaHyPar][Mt-KaHyPar] framework. There are currently no plans to merge this branch into the main Mt-KaHyPar project. If you plan to partition non-dynamic hypergraphs, we recommend using Mt-KaHyPar instead of DynPart.

### Reproducing Results
The results presented in the paper "DynPart: Fully Dynamic Hypergraph Partitioning" can be reproduced using the instructions below. The datasets and the raw benchmark result can be found in the [Zenodo repository](TODO).

Getting Started
-----------

Requirements:

 - A 64-bit Linux operating system.
 - A modern, C++17-ready compiler such as `g++` version 7 or higher, `clang` version 11.0.3 or higher, or `MinGW` compiler on Windows.
 - The [cmake][cmake] build system (>= 3.26).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).
   If you don't want to install boost by yourself, you can add the `-DKAHYPAR_DOWNLOAD_BOOST=On` flag
   to the cmake command to download, extract, and build the necessary dependencies automatically.
 - The [Intel Thread Building Blocks][tbb] library (TBB, minimum required version is OneTBB 2021.5.0).
   If you don't want to install TBB by yourself, you can add the `-DKAHYPAR_DOWNLOAD_TBB=On` flag
   to the cmake command to download oneTBB and extract the necessary dependencies automatically.
 - The [Portable Hardware Locality][hwloc] library (hwloc). This dependency is not strictly necessary,
   you can add the `-DKAHYPAR_DISABLE_HWLOC=On` flag to remove it.

The following command will install most of the required dependencies on a Ubuntu machine:

    sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev

### Build Commands

To build DynPart, use the following commands:

1. Clone the repository including submodules:

   ```git clone https://github.com/kahypar/dynpart.git```
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. --preset=default`
4. Run make: `make MtKaHyPar -j`

The build produces the executable `MtKaHyPar`, which can be found in `build/mt-kahypar/application/`.

Input Format
-----------

### Input Hypergraph Format
The algorithm expects an input hypergraph in the **hMetis** format. The first line of the input file must contain two integers: the number of vertices and the number of hyperedges. Each subsequent line represents a hyperedge and contains a list of vertex IDs that are part of that hyperedge. IDs are 1-based.

### Input Changelist Format
The algorithm expects a changelist file in the following format:
- The first line contains a single integer `n`, which is the **number of changes** in the changelist.
- Each Change is represented by 6 lines:
  - **Added Vertices**: Space separated list of ID's for vertices added to the hypergraph in the current change
  - **Added Hyperedges**: Space separated list of ID's for hyperedges added to the hypergraph
  - **Added Pins**: Space separated list of pairs of vertex and hyperedge ID's for pins added to the hypergraph, e.g. `(vertexID,hyperedgeID) (vertexID,hyperedgeID) ...`
  - **Removed Vertices**: Space separated list of ID's
  - **Removed Hyperedges**: Space separated list of ID's for hyperedges removed from the hypergraph
  - **Removed Pins**: Space separated list of pairs of vertex and hyperedge ID's for pins removed from the hypergraph, similar to the added pins format.

Running DynPart
-----------

To dynamically partition a **hypergraph** and its **changelist**, run the following command:

    ./mt-kahypar/application/MtKaHyPar -h <path-to-hgr> --d-changes-file <path-to-changelist-file> --preset-type=mutable -o km1 -k <# blocks> -e <imbalance (e.g. 0.03)> -v 0 --p-enable-community-detection=0 --d-strategy "localFM_rebalance_vcycle" --d-sb-threshold 5 --d-fm-buffer 100 --d-vcycle-step-size-pct <frequency of vcycles (e.g. 0.1)>

Run `--help` for a detailed description of the different program options. The most relevant options for dynamic partitioning are listed under **Dynamic Options** and start with the prefix `--d-`.

[//]: # (To enable writing the partition to a file after partitioning, you can add the following command line parameters to the partitioning call:)

[//]: # ()
[//]: # (    --d-write-partition-file=true --d-partition-output-folder=<path/to/folder>)

[//]: # ()
[//]: # (The partition file name is generated automatically based on parameters such as `k`, `imbalance`, `seed` and the input file name and will be located in the folder specified by `--partition-output-folder`. If you do not provide a partition output folder, the partition file will be placed in the same folder as the input hypergraph file.)

Reproducing Results
-----------

In the paper "DynPart: Fully Dynamic Hypergraph Partitioning" we benchmark DynPart in two different types of experiments. The first category evaluates the performance of DynPart on a set of 100 hypergraph partitioning benchmark instances with
three different fully dynamic update workloads. The second category evaluates the performance of DynPart on a set of graphs commonly used in dynamic graph partitioning, streaming the vertices.

### Fully Dynamic Hypergraph Partitioning Experiments

To reproduce the results of the fully dynamic hypergraph partitioning experiments, download the three changelist datasets from the [Zenodo repository](TODO). Each dataset contains the pre-processed hypergraphs and the corresponding changelists. More information on the individual changelist types can be found in the paper. The files corresponding to the datasets are:

1. **Mixed**: dynamic_hypergraph_changelist_mixed.zip
2. **Random Walk**: dynamic_hypergraph_changelist_random_walk.zip
3. **RTA**: dynamic_hypergraph_changelist_rta.zip

(Each of the experiments listed below has its dedicated result csv file in the Zenodo repository.)

All experiments sweep over the following parameters:
- -h: [path to hypergraph file]
- --d-changes-file: [path to changelist file corresponding to the hypergraph]
- -seed: [0, 42, 31415, 271828, 1618033]
- -k: [4,10,32,128]

Additionally the experiments, sweep over the following strategy configurations:

**Main Experiment**
- Greedy
  - --d-strategy: "greedy"
- Greedy + V-cycle 10% 
  - --d-strategy: "greedy_vcycle"
  - --d-vcycle-step-size-pct 0.1
- Greedy + V-cycle 5%
    - --d-strategy: "greedy_vcycle"
    - --d-vcycle-step-size-pct 0.05
- Static Repartitioning
    - --d-strategy: "repartition"
- Incremental Refiner
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-stop-vcycle: 0
    - --d-fm-buffer: -1
- Localized FM
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-stop-vcycle: 0
- Localized FM + V-cycle 10%
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-vcycle-step-size-pct 0.1
- Localized FM + V-cycle 5%
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-vcycle-step-size-pct 0.05

**Parameter Sweep Experiments**
- **Small Blocks Threshold**: Sweeps over --d-sb-threshold [0, 1, 2, 5, 10, 20, 100, 2147483647]
- **FM Buffer Size**: Sweeps over --d-fm-buffer [0, 5, 10, 100, 1000, 10000]
- **V-Cycle Stepsize**: Sweeps over --d-vcycle-step-size-pct [0.01, 0.05, 0.1, 0.2]

**Other Experiments**

The other reported experiments like the analysis of sub-runtimes are performed using the data from the main experiment.

### Streaming Experiments

To reproduce the results of the streaming experiments, download the streaming graph datasets (streaming.zip) from the [Zenodo repository](TODO). Each dataset contains the pre-processed graph with the first 100 vertices added and the corresponding changelist streaming the remaining vertices.

The benchmarked results can be found in the Zenodo repository (streaming.csv).

The competitors [Mt-KaHyPar][Mt-KaHyPar] and [FREIGHT][FREIGHT] were benchmarked in their original implementations. Please refer to their respective repositories for instructions on how to build and run them. Two other competitors, [Leopard][Leopard] and [Hermes][Hermes], are implemented in the DynPart framework and can be run using the same command line interface as DynPart.

The following settings were used for the streaming experiments:
- -h: [path to graph file]
- --d-changes-file: [path to changelist file corresponding to the graph]
- -seed: [0, 42, 31415, 271828, 1618033]
- -k: [10]
- --stream-changes: 0 // This pre-loads all changes into memory before partitioning. This may effect the runtime for the faster strategies like greedy. Setting this option to 1 will stream the changes from disk and may increase the runtime but can be necessary for instances with a large memory footprint.

**Configs**

- Greedy
    - --d-strategy: "greedy"
- Incremental Refiner
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-stop-vcycle: 0
    - --d-fm-buffer: -1
- Localized FM
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-stop-vcycle: 0
- Localized FM + V-cycle 100%
    - --d-strategy: "localFM_rebalance_vcycle"
    - --d-vcycle-step-size-pct 1
- Hermes 10%
    - --d-strategy: "hermes"
    - --d-vcycle-step-size-pct 0.1
- Leopard
    - --d-strategy: "leopard"

Licensing
---------

DynPart is a free software provided under the MIT License.
For more information see the [LICENSE file][LF].
We distribute this framework freely to foster the use and development of hypergraph partitioning tools.
If you use Mt-KaHyPar in an academic setting please cite the appropriate papers.

    // Mt-KaHyPar-D
    @inproceedings{MT-KAHYPAR-D,
      title     = {Scalable Shared-Memory Hypergraph Partitioning},
      author    = {Gottesbüren, Lars and
                   Heuer, Tobias and
                   Sanders, Peter and
                   Schlag, Sebastian},
      booktitle = {23rd Workshop on Algorithm Engineering and Experiments (ALENEX 2021)},
      pages     = {16--30},
      year      = {2021},
      publisher = {SIAM},
      doi       = {10.1137/1.9781611976472.2},
    }

    // Mt-KaHyPar-Q
    @inproceedings{MT-KAHYPAR-Q,
      title     = {Shared-Memory $n$-level Hypergraph Partitioning},
      author    = {Lars Gottesb{\"{u}}ren and
                   Tobias Heuer and
                   Peter Sanders and
                   Sebastian Schlag},
      booktitle = {24th Workshop on Algorithm Engineering and Experiments (ALENEX 2022)},
      year      = {2022},
      publisher = {SIAM},
      month     = {01},
      doi       = {10.1137/1.9781611977042.11}
    }

    // Mt-KaHyPar-Q-F
    @inproceedings{MT-KaHyPar-Q-F,
      title       =	{Parallel Flow-Based Hypergraph Partitioning},
      author      =	{Lars Gottesb\"{u}ren and
                     Tobias Heuer and
                     Peter Sanders},
      booktitle   =	{20th International Symposium on Experimental Algorithms (SEA 2022)},
      pages       =	{5:1--5:21},
      year        =	{2022},
      volume      =	{233},
      publisher   =	{Schloss Dagstuhl -- Leibniz-Zentrum f{\"u}r Informatik},
      doi         =	{10.4230/LIPIcs.SEA.2022.5}
    }

    // Deterministic Partitioning
    @inproceedings{MT-KAHYPAR-SDET,
      author    = {Lars Gottesb{\"{u}}ren and
                   Michael Hamann},
      title     = {Deterministic Parallel Hypergraph Partitioning},
      booktitle = {European Conference on Parallel Processing (Euro-Par)},
      volume    = {13440},
      pages     = {301--316},
      publisher = {Springer},
      year      = {2022},
      doi       = {10.1007/978-3-031-12597-3\_19},
    }

    // Unconstrained Refinement
    @inproceedings{MT-KAHYPAR-UNCONSTRAINED,
      author       = {Nikolai Maas and
                      Lars Gottesb{\"{u}}ren and
                      Daniel Seemaier},
      editor       = {Rezaul Chowdhury and
                      Solon P. Pissis},
      title        = {Parallel Unconstrained Local Search for Partitioning Irregular Graphs},
      booktitle    = {Symposium on Algorithm Engineering and Experiments (ALENEX 2024)},
      pages        = {32--45},
      publisher    = {{SIAM}},
      year         = {2024},
      doi          = {10.1137/1.9781611977929.3},
    }

    // Steiner Tree Objective
    @inproceedings{MT-KAHYPAR-STEINER-TREES,
      author       = {Tobias Heuer},
      editor       = {Rezaul Chowdhury and
                      Solon P. Pissis},
      title        = {A Direct \emph{k-}Way Hypergraph Partitioning Algorithm for Optimizing
                      the Steiner Tree Metric},
      booktitle    = {Symposium on Algorithm Engineering and Experiments (ALENEX 2024)},
      pages        = {15--31},
      publisher    = {{SIAM}},
      year         = {2024},
      doi          = {10.1137/1.9781611977929.2}
    }

    // Dissertation of Lars Gottesbüren
    @phdthesis{MT-KAHYPAR-DIS-GOTTESBUEREN,
      author         = {Lars Gottesb\"{u}ren},
      year           = {2023},
      title          = {Parallel and Flow-Based High-Quality Hypergraph Partitioning},
      doi            = {10.5445/IR/1000157894},
      pagetotal      = {256},
      school         = {Karlsruhe Institute of Technology}
    }

    // Dissertation of Tobias Heuer
    @phdthesis{MT-KAHYPAR-DIS-HEUER,
        author       = {Heuer, Tobias},
        year         = {2022},
        title        = {Scalable High-Quality Graph and Hypergraph Partitioning},
        doi          = {10.5445/IR/1000152872},
        pagetotal    = {242},
        school       = {Karlsruhe Institute of Technology}
    }

    // Mt-KaHyPar Journal Paper
    @article{MT-KAHYPAR-JOURNAL,
      author       = {Lars Gottesb{\"{u}}ren and
                      Tobias Heuer and
                      Nikolai Maas and
                      Peter Sanders and
                      Sebastian Schlag},
      title        = {Scalable High-Quality Hypergraph Partitioning},
      journal      = {{ACM} Transactions on Algorithms},
      volume       = {20},
      number       = {1},
      pages        = {9:1--9:54},
      year         = {2024},
      doi          = {10.1145/3626527},
    }

[cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
[tbb]: https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html
[hwloc]: https://www.open-mpi.org/projects/hwloc/
[Mt-KaHyPar]: https://github.com/kahypar/mt-kahypar
[FREIGHT]: https://github.com/KaHIP/FREIGHT
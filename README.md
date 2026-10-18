<h1 align="center">Dyn-HyPar - Dynamic Graph and Hypergraph Partitioner</h1>

License|Zenodo
:--:|:--:
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)|[![DOI](https://zenodo.org/badge/205879380.svg)](https://zenodo.org/badge/latestdoi/205879380)

Table of Contents
-----------

   * [Dyn-HyPar](#dyn-hypar)
   * [Getting Started](#getting-started)
   * [Build Commands](#build-commands)
   * [Running Dyn-HyPar](#running-dyn-hypar)
   * [Write Partition to Output File](#write-partition-to-output-file)
   * [Licensing](#licensing)


Dyn-HyPar
-----------
Dyn-HyPar is a dynamic hypergraph partitioning algorithm that can efficiently maintain a high-quality balanced partition of a hypergraph under dynamic changes. Changes include the insertion and deletion of hyperedges, vertices and pins. Dyn-HyPar computes a initial partition of the hypergraph using the Mt-KaHyPar framework. It utilizes incremental refinement, localized FM-based refinement and V-Cycles to move vertices between blocks on the fly and maintain a high-quality partition.

### Mt-KaHyPar
This branch features the dynamic hypergraph partitioning algorithm Dyn-HyPar, which is based on the Mt-KaHyPar framework. There are currently no plans to merge this branch into the main Mt-KaHyPar branch. If you plan to partition non-dynamic hypergraphs, we recommend using the main branch instead.

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

To build Dyn-HyPar, use the following commands:

1. Clone the repository including submodules:

   ```git clone https://github.com/kahypar/mt-kahypar.git```
2. Switch to the `dyn-hypar` branch: `git checkout dynamic-datastructure`
3. Create a build directory: `mkdir build && cd build`
5. Run cmake: `cmake .. --preset=default`
6. Run make: `make MtKaHyPar -j`

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

Running Dyn-HyPar
-----------

To dynamically partition a **hypergraph** and its **changelist**, run the following command:

    ./mt-kahypar/application/MtKaHyPar -h <path-to-hgr> --d-changes-file <path-to-changelist-file> --preset-type=mutable -o km1 -k <# blocks> -e <imbalance (e.g. 0.03)> -v 0 --p-enable-community-detection=0 --d-strategy "localFM_rebalance_vcycle" --d-sb-threshold 5 --d-fm-buffer 100 --d-vcycle-step-size-pct <frequency of vcycles (e.g. 0.1)>

Run `--help` for a detailed description of the different program options. The most relevant options for dynamic partitioning are listed under **Dynamic Options** and start with the prefix `--d-`.

### Write Partition to Output File

To enable writing the partition to a file after partitioning, you can add the following command line parameters to the partitioning call:

    --d-write-partition-file=true --d-partition-output-folder=<path/to/folder>

The partition file name is generated automatically based on parameters such as `k`, `imbalance`, `seed` and the input file name and will be located in the folder specified by `--partition-output-folder`. If you do not provide a partition output folder, the partition file will be placed in the same folder as the input hypergraph file.

Licensing
---------

Mt-KaHyPar is a free software provided under the MIT License.
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
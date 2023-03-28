<h1 align="center">Mt-KaHyPar - Multi-Threaded Karlsruhe Hypergraph Partitioner</h1>

License|Linux & Windows Build|Code Coverage|Zenodo
:--:|:--:|:--:|:--:
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)|[![Build Status](https://github.com/kahypar/mt-kahypar/actions/workflows/mt_kahypar_ci.yml/badge.svg)](https://github.com/kahypar/mt-kahypar/actions/workflows/mt_kahypar_ci.yml)|[![codecov](https://codecov.io/gh/kahypar/mt-kahypar/branch/master/graph/badge.svg?token=sNWRRtXZjI)](https://codecov.io/gh/kahypar/mt-kahypar)|[![DOI](https://zenodo.org/badge/205879380.svg)](https://zenodo.org/badge/latestdoi/205879380)

Table of Contents
-----------

   * [About Mt-KaHyPar](#about-mt-kahypar)
   * [Requirements](#requirements)
   * [Building Mt-KaHyPar](#building-mt-kahypar)
   * [Running Mt-KaHyPar](#running-mt-kahypar)
   * [Performance](#performance)
   * [The C Library Interface](#the-c-library-interface)
   * [The Python Library Interface](#the-python-library-interface)
   * [Licensing](#licensing)

About Mt-KaHyPar
-----------
Mt-KaHyPar is a shared-memory algorithm for partitioning graphs and hypergraphs. The balanced (hyper)graph partitioning problem
asks for a partition of the node set of a (hyper)graph into *k* disjoint blocks of roughly the same size (usually a small imbalance
is allowed by at most 1 + ε times the average block weight), while simultanously minimizing an objective function defined on the (hyper)edges.
The edge-cut metric is the most prominent objective function for graph partitioning, which sums over the weight of all edges that connect
two blocks. For hypergraph partitioning, research has focused on the connectivity metric that additionally multiplies the weight of each
hyperedge with the number of blocks connected by that hyperedge (sum over the terms (λ(e) − 1) * ω(e) where λ(e) is the number of blocks connected by hyperedge e and ω(e) is the weight of the hyperedge).

When we started to work on this topic, we realized there was a large gap between the solution quality of the partitions produced by sequential and parallel partitioning algorithms. We then started to parallelize all techniques used in the best sequential partitioning algorithms without compromises in solution quality. The main outcome of our work is a parallel partitioning algorithm that can partition extremely large graphs and hypergraphs (with billion of edges) with comparable solution quality to the best sequential graph partitioner [KaFFPa](https://github.com/KaHIP/KaHIP) and hypergraph partitioner [KaHyPar](https://kahypar.org/) while being (more) than an order of magnitude faster with only ten threads.

Initially, we focused on hypergraph partitioning but recently implemented optimized data structures for graph partitioning (which led to a speedup by a factor of two for graphs). Besides our high-quality configuration, we provide several other faster configurations that are already
able to outperform most of the existing partitioning algorithms with regard to solution quality and running time. Moreover, we also provide a deterministic version of our partitioning algorithm. We refer the reader to our [publications](#licensing) for more information.

<img src="https://cloud.githubusercontent.com/assets/484403/25314222/3a3bdbda-2840-11e7-9961-3bbc59b59177.png" alt="alt text" width="50%" height="50%"><img src="https://cloud.githubusercontent.com/assets/484403/25314225/3e061e42-2840-11e7-860c-028a345d1641.png" alt="alt text" width="50%" height="50%">

Requirements
-----------

The Multi-Threaded Karlsruhe Hypergraph Partitioning Framework requires:

  - A 64-bit Linux or Windows operating system.
  - A modern, ![C++17](https://img.shields.io/badge/C++-17-blue.svg?style=flat)-ready compiler such as `g++` version 7 or higher, `clang` version 11.0.3 or higher, or `MinGW` compiler on Windows (tested with version 12.1).
 - The [cmake][cmake] build system (>= 3.16).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).
   If you don't want to install boost by yourself, you can add the `-DKAHYPAR_DOWNLOAD_BOOST=On` flag
   to the cmake command to download, extract, and build the neccessary dependencies automatically.
 - The [Intel Thread Building Blocks][tbb] library (TBB, minimum required version is OneTBB 2021.5.0).
   If you don't want to install TBB by yourself, you can add the `-DKAHYPAR_DOWNLOAD_TBB=On` flag (only available on Linux)
   to the cmake command to download oneTBB 2021.7.0 and extract the neccessary dependencies automatically.
   Mt-KaHyPar also compiles with older version of TBB. However, we observed unexpected behavior of a TBB function
   on which we rely on which causes on our side a segmentation fault in really rare cases. If you want to ignore these
   warning, you can add `-DKAHYPAR_ENFORCE_MINIMUM_TBB_VERSION=OFF` to the cmake build command.
 - The [Portable Hardware Locality][hwloc] library (hwloc)

### Linux

The following command will install most of the required dependencies on a Ubuntu machine:

    sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev

### Windows

The following instructions setup the environment used to build Mt-KaHyPar on Windows machines:

  1. Download and install [MSYS2][MSYS2] from the official website (https://www.msys2.org/).
  2. Launch the `MSYS2 MinGW x64` terminal.
  3. Update the package manager database by running the following command:

    pacman -Syu

  4. The following command will then install all required dependencies:

    pacman -S make mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc mingw-w64-x86_64-python3 mingw-w64-x86_64-tbb

  5. Rename `libtbb12.dll.a` to `libtbb.dll.a` which is located in `C:\msys64\mingw64\lib` (or `/mingw64/lib` within the
     `MSYS2 MinGW x64` terminal)

Please **note** that Mt-KaHyPar was primarily tested and evaluated on Linux machines. While a Windows build has been provided and tested on `MSYS2` using `pacman` to install the required dependencies, we cannot provide any performance guarantees or ensure that the Windows version is free of bugs. At this stage, Windows support is experimental. We are happy to accept contributions to improve Windows support.

Building Mt-KaHyPar
-----------

To build Mt-KaHyPar, you can run the `build.sh` script (creates a `build` folder) or use the following commands:

1. Clone the repository including submodules:

   ```git clone --depth=1 --recursive git@github.com:kahypar/mt-kahypar.git```

2. Create a build directory: `mkdir build && cd build`
3. *Only on Windows machines*: `export CMAKE_GENERATOR="MSYS Makefiles"`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE` (on Windows machines add `-DKAHYPAR_DOWNLOAD_BOOST=On`)
4. Run make: `make MtKaHyPar -j`

The build produces the executable `MtKaHyPar`, which can be found in `build/mt-kahypar/application/`.

Running Mt-KaHyPar
-----------

Mt-KaHyPar has several configuration parameters. We recommend to use one of our presets (also located in the `config` folder):

- `default`: corresponds to Mt-KaHyPar-D in our publications (`config/default_preset.ini`)
- `default_flows`: corresponds to Mt-KaHyPar-D-F (`config/default_flow_preset.ini`)
- `quality`: corresponds to Mt-KaHyPar-Q (`config/quality_preset.ini`)
- `quality_flows`: corresponds to Mt-KaHyPar-Q-F (`config/quality_flow_preset.ini`)
- `large_k`: configuration for partitioning a (hyper)graph into a large number of blocks (e.g. >= 1024 blocks, `config/large_k_preset.ini`)
- `deterministic`: configuration for deterministic partitioning (`config/deterministic_preset.ini`)

The presets can be ranked from lowest to the highest quality as follows: `large_k`, `deterministic`,
`default`, `quality`, `default_flows` and `quality_flows`.
Initially, we started with the `default` and `quality` configuration and then extended both configurations with flow-based refinement (`default_flows` and `quality_flows`). We then found that our `default_flows` configuration produces better partitions than the `quality` configuration. However, we still keep the naming due to the naming in our publications. In general, we recommend to use the `default` configuration to compute good partitions very fast and the `default_flows` configuration to compute high-quality solutions. The `quality_flows` configuration computes better partitions than our `default_flows` configuration by 0.5% on average at the cost of a two times longer running time for medium-sized instances (up to 100 million pins). When you have to partition a (hyper)graph into a large number of blocks (e.g., >= 1024 blocks), you can use our `large_k` configuration. However, we only recommend to use this if you experience high running times with one of our other configurations as this can significantly worsen the partitioning quality.

If you want to change configuration parameters manually, please run `--help` for a detailed description of the different program options. We use the [hMetis format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf) for hypergraph files as well as the partition output file and the [Metis format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf) for graph files. Per default, we expect the input to be in hMetis format, but you can read graphs in Metis format via command line parameter `--input-file-format=metis`. If your input file is a graph, you can switch to our optimized graph data structures via command line parameter `--instance-type=graph`.

To run Mt-KaHyPar, you can use the following command:

    ./MtKaHyPar -h <path-to-hgr> --preset-type=<large_k/deterministic/default/default_flows/quality/quality_flows> --instance_type=<hypergraph/graph> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1

or directly provide a configuration file (see `config` folder):

    ./MtKaHyPar -h <path-to-hgr> -p <path-to-config-file> --instance_type=<hypergraph/graph> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1

To enable writing the partition to a file set the flag `--write-partition-file=true`.
By default the file will be placed in the same folder as the input hypergraph file. Set `--partition-output-folder=path/to/folder` to specify a desired output folder. The partition file name is generated automatically based on parameters such as `k`, `imbalance`, `seed` and the input file name.

Further, there are several useful options that can provide you with additional insights during and after the partitioning process:
- `--verbose=true`: Displays detailed information on the partitioning process
- `--show-detailed-timings=true`: Shows detailed subtimings of each phase of the algorithm at the end of partitioning
- `--enable-progress-bar=true`: Shows a progess bar during the coarsening and refinement phase

Mt-KaHyPar uses 32-bit node and hyperedge IDs. If you want to partition hypergraphs with more than 4.294.967.295 nodes or hyperedges, add option `-DKAHYPAR_USE_64_BIT_IDS=ON` to the `cmake` build command.

Performance
-----------

We have summarized our experimental results on an [external webpage][ExperimentalResults]. The resource provides a detailed
overview of Mt-KaHyPar's performance compared to other prominent state-of-the-art systems in terms of running time
and quality.


The C Library Interface
-----------

We provide a simple C-style interface to use Mt-KaHyPar as a library.  The library can be built and installed via

```sh
make install.mtkahypar # use sudo (Linux) or run shell as an adminstrator (Windows) to install system-wide
```

Note: When installing locally, the build will exit with an error due to missing permissions.
However, the library is still built successfully and is available in the build folder.

The library interface can be found in `include/libmtkahypar.h` with a detailed documentation of its functionality. We also provide several examples in the folder `lib/examples` that show how to use the library.

Here is a short example how you can partition a hypergraph using our library interface:

```cpp
#include <memory>
#include <vector>
#include <iostream>
#include <thread>

#include <libmtkahypar.h>

int main(int argc, char* argv[]) {

  // Initialize thread pool
  mt_kahypar_initialize_thread_pool(
    std::thread::hardware_concurrency() /* use all available cores */,
    true /* activate interleaved NUMA allocation policy */ );

  // Setup partitioning context
  mt_kahypar_context_t* context = mt_kahypar_context_new();
  mt_kahypar_load_preset(context, DEFAULT /* corresponds to MT-KaHyPar-D */);
  // In the following, we partition a hypergraph into two blocks
  // with an allowed imbalance of 3% and optimize the connective metric (KM1)
  mt_kahypar_set_partitioning_parameters(context,
    2 /* number of blocks */, 0.03 /* imbalance parameter */,
    KM1 /* objective function */, 42 /* seed */);
  // Enable logging
  mt_kahypar_set_context_parameter(context, VERBOSE, "1");

  // Load Hypergraph for DEFAULT preset
  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar_read_hypergraph_from_file(
      "path/to/hypergraph/file", DEFAULT, HMETIS /* file format */);

  // Partition Hypergraph
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_partition(hypergraph, context);

  // Extract Partition
  std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
    std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
  mt_kahypar_get_partition(partitioned_hg, partition.get());

  // Extract Block Weights
  std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
    std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
  mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

  // Compute Metrics
  const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
  const double km1 = mt_kahypar_km1(partitioned_hg);

  // Output Results
  std::cout << "Partitioning Results:" << std::endl;
  std::cout << "Imbalance         = " << imbalance << std::endl;
  std::cout << "Km1               = " << km1 << std::endl;
  std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
  std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;

  mt_kahypar_free_context(context);
  mt_kahypar_free_hypergraph(hypergraph);
  mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}
```

To compile the program using `g++` run:

```sh
g++ -std=c++17 -DNDEBUG -O3 your_program.cc -o your_program -lmtkahypar
```

To execute the binary, you need to ensure that the installation directory
(probably `/usr/local/lib` (Linux) and `C:\Program Files (x86)\MtKaHyPar\bin` (Windows) for system-wide installation)
is included in the dynamic library path.
The path can be updated on Linux with:

```sh
LD_LIBRARY_PATH="$LD_LIBRARY_PATH;/usr/local/lib"
export LD_LIBRARY_PATH
```

On Windows, add `C:\Program Files (x86)\KaHyPar\bin` to `PATH` in the environment variables settings.

To remove the library from your system use the provided uninstall target:

```sh
make uninstall-mtkahypar
```

**Note** that we internally use different data structures to represent a (hyper)graph based on the corresponding configuration (`mt_kahypar_preset_type_t`). The `mt_kahypar_hypergraph_t` structure stores a pointer to this data structure and also a type description. Therefore, you can not partition a (hyper)graph with all available configurations once it is loaded or constructed. However, you can check the compatibility of a hypergraph with a configuration with the following code:

```cpp
mt_kahypar_context_t context = mt_kahypar_context_new();
mt_kahypar_load_preset(context, QUALITY);
// Check if the hypergraph is compatible with the QUALITY preset
if ( check_compatibility(hypergraph, QUALITY) ) {
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_partition(hypergraph, context);
}
```

The Python Library Interface
-----------

You can install the Python library interface via

```sh
make mtkahypar_python
```

This will create a shared library in the `build/python` folder (`mtkahypar.so` on Linux and `mtkahypar.pyd` on Windows).
Copy the libary to your Python project directory to import Mt-KaHyPar as a Python module.

A documentation of the Python module can be found in `python/module.cpp`, or by importing the module (`import mtkahypar`) and calling `help(mtkahypar)` in Python. We also provide several examples that show how to use the Python interface in the folder `python/examples`.

Here is a short example how you can partition a hypergraph using our python interface:

```py
import multiprocessing
import mtkahypar

# Initialize thread pool
mtkahypar.initializeThreadPool(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtkahypar.Context()
context.loadPreset(mtkahypar.PresetType.DEFAULT) # corresponds to Mt-KaHyPar-D
# In the following, we partition a hypergraph into two blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.setPartitioningParameters(
  2,                       # number of blocks
  0.03,                    # imbalance parameter
  mtkahypar.Objective.KM1, # objective function
  42)                      # seed
context.logging = True     # enables partitioning output

# Load hypergraph from file
hypergraph = mtkahypar.Hypergraph(
  "path/to/hypergraph/file", # hypergraph file
  mtkahypar.FileFormat.HMETIS) # hypergraph is stored in hMetis file format

# Partition hypergraph
partitioned_hg = hypergraph.partition(context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance()))
print("km1       = " + str(partitioned_hg.km1()))
print("Block Weights:")
print("Weight of Block 0 = " + str(partitioned_hg.blockWeight(0)))
print("Weight of Block 1 = " + str(partitioned_hg.blockWeight(1)))
```

We also provide an optimized graph data structure for partitioning plain graphs. The following example loads and partitions a graph:

```py
# Load graph from file
graph = mtkahypar.Graph(
  "path/to/graph/file", # graph file
  mtkahypar.FileFormat.METIS) # graph is stored in Metis file format

# Partition graph
partitioned_graph = graph.partition(context)
```
**Note** that when you want to partition a hypergraph into large number of blocks (e.g., k > 1024), you can use our `LARGE_K` confguration and the `partitionIntoLargeK(...)` function of the hypergraph object. If you use an other configuration for large k partitioning, you may run into memory and running time issues during partitioning. However, this depends on the size of the hypergraph and the memory capacity of your target machine. For partitioning plain graphs, you can load the `LARGE_K` configuration, but you can still use the `partition(...)` function of the graph object. Here is an example that partitions a hypergraph into 1024 blocks:

```py
# Setup partitioning context
context = mtkahypar.Context()
context.loadPreset(mtkahypar.PresetType.LARGE_K)
# In the following, we partition a hypergraph into 1024 blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.setPartitioningParameters(1024, 0.03, mtkahypar.Objective.KM1, 42)

# Load and partition hypergraph
hypergraph = mtkahypar.Hypergraph("path/to/hypergraph/file", mtkahypar.FileFormat.HMETIS)
partitioned_hg = hypergraph.partitionIntoLargeK(context)
```

Bug Reports
-----------

We encourage you to report any problems with Mt-KaHyPar via the [github issue tracking system](https://github.com/kittobi1992/mt-kahypar/issues) of the project.

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

    // Deterministic Partitioning (Technical Report)
    @techreport{MT-KAHYPAR-SDET,
      title       = {Deterministic Parallel Hypergraph Partitioning},
      author      = {Lars Gottesbüren and
                     Michael Hamann},
      institution = {Karlsruhe Institute of Technology},
      year        = {2021},
      url         = {https://arxiv.org/pdf/2112.12704.pdf}
    }

Contributing
------------
If you are interested in contributing to the Mt-KaHyPar framework
feel free to contact us or create an issue on the
[issue tracking system](https://github.com/kahypar/mt-kahypar/issues).

[cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
[tbb]: https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html
[hwloc]: https://www.open-mpi.org/projects/hwloc/
[LF]: https://github.com/kahypar/mt-kahypar/blob/master/LICENSE "License"
[SetA]: http://algo2.iti.kit.edu/heuer/alenex21/instances.html?benchmark=set_a
[SetB]: http://algo2.iti.kit.edu/heuer/alenex21/instances.html?benchmark=set_b
[ExperimentalResults]: https://algo2.iti.kit.edu/heuer/mt_kahypar/
[MSYS2]: https://www.msys2.org/

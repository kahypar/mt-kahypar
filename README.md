<h1 align="center">Mt-KaHyPar - Multi-Threaded Karlsruhe Graph and Hypergraph Partitioner</h1>

License|Linux, MacOS & Windows Build|Code Coverage|Zenodo
:--:|:--:|:--:|:--:
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)|[![Build Status](https://github.com/kahypar/mt-kahypar/actions/workflows/master_ci.yml/badge.svg)](https://github.com/kahypar/mt-kahypar/actions/workflows/mt_kahypar_ci.yml)|[![codecov](https://codecov.io/gh/kahypar/mt-kahypar/branch/master/graph/badge.svg?token=sNWRRtXZjI)](https://codecov.io/gh/kahypar/mt-kahypar)|[![DOI](https://zenodo.org/badge/205879380.svg)](https://zenodo.org/badge/latestdoi/205879380)

Table of Contents
-----------

   * [About Mt-KaHyPar](#about-mt-kahypar)
   * [Features](#features)
   * [Installing Mt-KaHyPar](#installing-mt-kahypar)
   * [Building Mt-KaHyPar from Source](#building-mt-kahypar-from-source)
   * [Running Mt-KaHyPar](#running-mt-kahypar)
   * [Using Mt-KaHyPar as a library](#using-mt-kahypar-as-a-library)
   * [Supported Objective Functions](#supported-objective-functions)
   * [Improving Compile Times](#improving-compile-times)
   * [Licensing](#licensing)

About Mt-KaHyPar
-----------
Mt-KaHyPar is a shared-memory algorithm for partitioning graphs and hypergraphs. The balanced (hyper)graph partitioning problem
asks for a partition of the node set of a (hyper)graph into *k* disjoint blocks of roughly the same size (usually a small imbalance
is allowed by at most 1 + ε times the average block weight), while simultaneously minimizing an objective function defined on the (hyper)edges.
Mt-KaHyPar can optimize the cut-net, connectivity, sum-of-external-degrees, and Steiner tree metric (see [Supported Objective Functions](#supported-objective-functions)).

<img src="https://cloud.githubusercontent.com/assets/484403/25314222/3a3bdbda-2840-11e7-9961-3bbc59b59177.png" alt="alt text" width="50%" height="50%"><img src="https://cloud.githubusercontent.com/assets/484403/25314225/3e061e42-2840-11e7-860c-028a345d1641.png" alt="alt text" width="50%" height="50%">

The highest-quality configuration of Mt-KaHyPar computes partitions that are on par with those produced by the best sequential partitioning algorithms, while being almost an order of magnitude faster with only *ten* threads (e.g., when compared to [KaFFPa](https://github.com/KaHIP/KaHIP) or [KaHyPar](https://kahypar.org/)).
Besides our high-quality configuration, we provide several other faster configurations that are already able to outperform most of the existing partitioning algorithms with regard to solution quality and running time.
The figure below summarizes the time-quality trade-off of different hypergraph (left, connectivity metric) and graph partitioning algorithms (right, cut-net metric).
The plot is based on an experiment with over 800 graphs and hypergraphs and relates the average solution quality and running time of each algorithm to the best achievable results.
Points on the lower-left are considered better. Partially transparent markers indicate solvers producing more than 15% infeasible partitions (either imbalanced or timeout).
For more details, we refer the reader to our [publications](#licensing).

![time_quality_trade_off](https://github.com/kahypar/mt-kahypar/assets/9654047/a5cc1c41-5ca5-496a-ba50-91965e73226b)

Features
-----------

Besides its fast and high-quality partitioning algorithm, Mt-KaHyPar provides many other useful features:

- **Scalability**: Mt-KaHyPar has excellent scaling behaviour (up to 25 with 64 threads), while increasing the number of threads does not adversely affect the solution quality.
- **Deterministic Partitioning**: Mt-KaHyPar offers a deterministic partitioning algorithm, ensuring consistent solutions for the same input and random seed.
- **Large K Partitioning**: We provide a partitioning configuration for partitioning (hyper)graphs into a large number of blocks (e.g., k > 1024).
- **Graph Partitioning**: Mt-KaHyPar includes optimized data structures for graph partitioning, achieving a speedup by a factor of two for plain graphs.
- **Objective Functions**: Mt-KaHyPar can optimize the cut-net, connectivity, and sum-of-external-degrees metric (for more details, see [Supported Objective Functions](#supported-objective-functions))
- **Mapping (Hyper)Graphs Onto Graphs**: In many applications (e.g., distributed computation), the partition is assigned to a communication network that can be represented as a graph.
  However, conventional objective functions do not consider the topology of the target graph.
  We therefore provide a mode that maps the nodes of a (hyper)graph onto a target graph via the [Steiner tree metric](#steiner-tree-metric).
- **Fixed Vertices**: Fixed vertices are nodes that are preassigned to a particular block and are not allowed to change their block during partitioning.


Installing Mt-KaHyPar
-----------

For Linux (x86) and MacOS, the [Mt-KaHyPar Python package](https://pypi.org/project/mtkahypar/) can be installed via pip:

    pip install mtkahypar


We also provide Debian packages that contain the CLI application and the C library interface for the [latest release](https://github.com/kahypar/mt-kahypar/releases/latest).


Building Mt-KaHyPar from Source
-----------

Mt-KaHyPar requires:

 - A 64-bit Linux, MacOS, or Windows operating system.
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

### Linux

The following command will install most of the required dependencies on a Ubuntu machine:

    sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev

### MacOS

The following command will install most of the required dependencies on a MacOS machine:

    brew install tbb boost hwloc

### Windows

The following instructions set up the environment used to build Mt-KaHyPar on Windows machines:

  1. Download and install [MSYS2][MSYS2] from the official website (https://www.msys2.org/).
  2. Launch the `MSYS2 MinGW x64` terminal.
  3. Update the package manager database by running the following command:

    pacman -Syu

  4. The following command will then install all required dependencies:

    pacman -S make mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc mingw-w64-x86_64-python3 mingw-w64-x86_64-boost mingw-w64-x86_64-tbb

### Build Commands

To build Mt-KaHyPar, use the following commands:

1. Clone the repository including submodules:

   ```git clone https://github.com/kahypar/mt-kahypar.git```

2. Create a build directory: `mkdir build && cd build`
3. *Only on Windows machines: `export CMAKE_GENERATOR="MSYS Makefiles"`*
3. Run cmake: `cmake .. --preset=<default/python/dev>`
4. Run make: `make MtKaHyPar -j`

The build produces the executable `MtKaHyPar`, which can be found in `build/mt-kahypar/application/`.

As a user of Mt-KaHyPar, the `default` cmake preset is appropriate (or `python` for installing the Python interface).
If you work on Mt-KaHyPar or want to run benchmarks, use the `dev` preset.

Please note that Mt-KaHyPar was primarily tested and evaluated on Linux machines. While a Windows build has been provided and tested on `MSYS2` using `pacman` to install the required dependencies, we cannot provide any performance guarantees or ensure that the Windows version is free of bugs. We are happy to accept contributions to improve Windows support.

Running Mt-KaHyPar
-----------

To partition a **hypergraph** with our default configuration, you can use the following command:

    ./mt-kahypar/application/MtKaHyPar -h <path-to-hgr> --preset-type=default -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1

### Partitioning Configurations

Mt-KaHyPar provides several partitioning configurations with different time-quality trade-offs. The configurations are stored in `ini` files located in the `config` folder. However, we recommend using the `--preset-type` command line parameter to run Mt-KaHyPar with a specific partitioning configuration:

    --preset-type=<large_k/deterministic/default/quality/highest_quality>

- `large_k`: configuration for partitioning (hyper)graphs into a large number of blocks (e.g. >= 1024 blocks)
- `deterministic`: configuration for deterministic partitioning
- `default`: computes good partitions very fast
- `quality`: computes high-quality partitions (uses flow-based refinement)
- `highest_quality`: highest-quality configuration (uses n-level coarsening and flow-based refinement)

The presets can be ranked from lowest to the highest-quality as follows: `large_k`, `deterministic`,
`default`, `quality`, and `highest_quality`.
We recommend using the `default` configuration to compute good partitions very fast and the `quality` configuration to compute high-quality solutions.
The `highest_quality` configuration computes better partitions than our `quality` configuration by 0.5% on average at the cost of a two times longer running time for medium-sized instances (up to 100 million pins).
When you have to partition a (hyper)graph into a large number of blocks (e.g., >= 1024 blocks), you can use our `large_k` configuration.
However, we only recommend using this if you experience high running times with one of our other configurations as this can significantly worsen the partitioning quality.

### Objective Functions

Mt-KaHyPar can optimize the cut-net, connectivity, and sum-of-external-degrees metric (see [Supported Objective Functions](#supported-objective-functions)).

    -o <cut/km1/soed>

### Graph Partitioning

To partition a **graph** with Mt-KaHyPar, you can add the following command line parameters to the partitioning call:

    -h <path-to-graph> --instance-type=graph --input-file-format=<metis/hmetis> -o cut

Mt-KaHyPar then uses optimized data structures for graph partitioning, which speeds up the partitioning time by a factor of two compared to our hypergraph partitioning code. Per default, we expect the input in [hMetis format](mt-kahypar/io/docs/FileFormats.md#hmetis-format-for-input-hypergraphs), but you can read graph files in [Metis format](mt-kahypar/io/docs/FileFormats.md#metis-format-for-input-graphs) via `--input-file-format=metis`.

### Mapping (Hyper)Graphs onto Graphs

To map a **(hyper)graph** onto a **target graph** with Mt-KaHyPar, you can add the following command line parameters to the partitioning call:

    -g <path-to-target-graph> -o steiner_tree

The target graph is expected to be in [Metis format](mt-kahypar/io/docs/FileFormats.md#metis-format-for-input-graphs). The nodes of the (hyper)graph are then mapped onto the nodes of the target graph, while optimizing the Steiner tree metric (see [Supported Objective Functions](#supported-objective-functions)).

### Fixed Vertices

Fixed vertices are nodes that are preassigned to particular block and are not allowed to change their block during partitioning. Mt-KaHyPar reads fixed vertices from a file in the [hMetis fix file format](mt-kahypar/io/docs/FileFormats.md#hmetis-fix-file-format), which can be provided via the following command line parameter:

    -f <path-to-fixed-vertex-file>

Note that fixed vertices are only supported in our `default`, `quality`, and `highest_quality` configurations.

### Individual Target Block Weights

Per default, Mt-KaHyPar enforces that the weight of each block must be smaller than the average block weight (weight of the hypergraph divided by the number of blocks) times (1 + ε). However, you can provide individual target block weights for each block via

    --part-weights=weight_of_block_0 weight_of_block_1 ... weight_of_block_k

Note that the sum of all individual target block weights must be larger than the total weight of all nodes.

### Write Partition to Output File

To enable writing the partition to a file after partitioning, you can add the following command line parameters to the partitioning call:

    --write-partition-file=true --partition-output-folder=<path/to/folder>

The partition file name is generated automatically based on parameters such as `k`, `imbalance`, `seed` and the input file name and will be located in the folder specified by `--partition-output-folder`. If you do not provide a partition output folder, the partition file will be placed in the same folder as the input hypergraph file.

### Other Useful Program Options

There are several useful options that can provide you with additional insights during and after the partitioning process:

- `--verbose=true`: Displays detailed information on the partitioning process
- `--show-detailed-timings=true`: Shows detailed sub-timings of each phase of the algorithm at the end of partitioning
- `--enable-progress-bar=true`: Shows a progress bar during the coarsening and refinement phase


If you want to change other configuration parameters manually, please run `--help` for a detailed description of the different program options.

Using Mt-KaHyPar as a library
-----------

We provide a simple C interface to use Mt-KaHyPar as a library, as well as a Python interface.
On Linux or MacOS, the C library can be built and installed via

```sh
make install-mtkahypar  # use sudo (Linux & MacOS) or run shell as an administrator (Windows) to install system-wide
```

Note: When installing locally, the build will exit with an error due to missing permissions.
However, the library is still built successfully and is available in the build folder.

To remove the library from your system use the provided uninstall target:

```sh
make uninstall-mtkahypar
```

### Integration via Cmake

If possible, the best way to integrate the C library is directly via cmake using the `MtKaHyPar::mtkahypar` target.

If the library is installed on the system, it can be used via `find_package`:

```cmake
find_package(MtKaHyPar)
if(MtKaHyPar_FOUND)
  add_executable(example example.cc)
  target_link_libraries(example MtKaHyPar::mtkahypar)
endif()
```

Alternatively, you can use Mt-KaHyPar directly via `FetchContent`:

```cmake
FetchContent_Declare(
  MtKaHyPar EXCLUDE_FROM_ALL
  GIT_REPOSITORY https://github.com/kahypar/mt-kahypar
  GIT_TAG        v1.5
)
FetchContent_MakeAvailable(MtKaHyPar)

add_executable(example example.cc)
target_link_libraries(example MtKaHyPar::mtkahypar)
```

When including Mt-KaHyPar directly, it is also possible to control static versus dynamic linking with the `BUILD_SHARED_LIBS` and `KAHYPAR_STATIC_LINK_DEPENDENCIES` cmake options
(note that static linking support is still experimental and not available for the installed library).

### The C Library Interface

The library interface can be found in [`include/mtkahypar.h`](include/mtkahypar.h) with a detailed documentation. We also provide [several examples](lib/examples) that show how to use the library.

Here is a short example of how you can partition a hypergraph using our library interface:

```cpp
#include <cassert>
#include <memory>
#include <vector>
#include <iostream>
#include <thread>

#include <mtkahypar.h>

int main(int argc, char* argv[]) {
  mt_kahypar_error_t error{};

  // Initialize
  mt_kahypar_initialize(
    std::thread::hardware_concurrency() /* use all available cores */,
    true /* activate interleaved NUMA allocation policy */ );

  // Setup partitioning context
  mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
  // In the following, we partition a hypergraph into two blocks
  // with an allowed imbalance of 3% and optimize the connective metric (KM1)
  mt_kahypar_set_partitioning_parameters(context,
    2 /* number of blocks */, 0.03 /* imbalance parameter */,
    KM1 /* objective function */);
  mt_kahypar_set_seed(42 /* seed */);
  // Enable logging
  mt_kahypar_status_t status =
    mt_kahypar_set_context_parameter(context, VERBOSE, "1", &error);
  assert(status == SUCCESS);

  // Load Hypergraph for DEFAULT preset
  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar_read_hypergraph_from_file("path/to/hypergraph/file",
      context, HMETIS /* file format */, &error);
  if (hypergraph.hypergraph == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }

  // Partition Hypergraph
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_partition(hypergraph, context, &error);
  if (partitioned_hg.partitioned_hg == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }

  // Extract Partition
  auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(
    mt_kahypar_num_hypernodes(hypergraph));
  mt_kahypar_get_partition(partitioned_hg, partition.get());

  // Extract Block Weights
  auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
  mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

  // Compute Metrics
  const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
  const int km1 = mt_kahypar_km1(partitioned_hg);

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

We recommend [integrating the library via cmake](#integration-via-cmake) into your project, to manage compiling and linking automatically.

However, it is also possible to directly compile the program using `g++`:

```sh
g++ -std=c++17 -DNDEBUG -O3 your_program.cc -o your_program -lmtkahypar
```

To execute the produced binary, you need to ensure that the installation directory
(probably `/usr/local/lib` on Linux and `C:\Program Files (x86)\MtKaHyPar\bin` on Windows)
is included in the dynamic library path (`LD_LIBRARY_PATH` on Linux, `PATH` on Windows).

**Note** that we internally use different data structures to represent a (hyper)graph based on the corresponding configuration (`mt_kahypar_preset_type_t`).
The `mt_kahypar_hypergraph_t` structure stores a pointer to this data structure and also a type description.
Therefore, you can not partition a (hyper)graph with all available configurations once it is loaded or constructed. However, you can check the compatibility of a hypergraph with a configuration with the following code:

```cpp
mt_kahypar_context_t* context = mt_kahypar_context_from_preset(QUALITY);
// Check if the hypergraph is compatible with the QUALITY preset
if ( mt_kahypar_check_compatibility(hypergraph, QUALITY) ) {
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_partition(hypergraph, context, &error);
}
```

### The Python Library Interface

You can install the Python library interface via

```sh
make mtkahypar_python
```

This will create a shared library in the `build/python` folder (`mtkahypar.so` on Linux and `mtkahypar.pyd` on Windows).
Copy the libary to your Python project directory to import Mt-KaHyPar as a Python module.

A documentation of the Python module can be found by importing the module (`import mtkahypar`) and calling `help(mtkahypar)` in Python.
We also provide [several examples](python/examples) that show how to use the Python interface.

Here is a short example of how you can partition a hypergraph using our Python interface:

```py
import multiprocessing
import mtkahypar

# Initialize
mtk = mtkahypar.initialize(multiprocessing.cpu_count()) # use all available cores

# Setup partitioning context
context = mtk.context_from_preset(mtkahypar.PresetType.DEFAULT)
# In the following, we partition a hypergraph into two blocks
# with an allowed imbalance of 3% and optimize the connectivity metric
context.set_partitioning_parameters(
  2,                       # number of blocks
  0.03,                    # imbalance parameter
  mtkahypar.Objective.KM1) # objective function
mtkahypar.set_seed(42)     # seed
context.logging = True

# Load hypergraph from file (assumes hMetis file format per default)
hypergraph = mtk.hypergraph_from_file("path/to/hypergraph/file", context)

# Partition hypergraph
partitioned_hg = hypergraph.partition(context)

# Output metrics
print("Partition Stats:")
print("Imbalance = " + str(partitioned_hg.imbalance(context)))
print("km1       = " + str(partitioned_hg.km1()))
print("Block Weights:")
for i in partitioned_hg.blocks():
  print(f"Weight of Block {i} = {partitioned_hg.block_weight(i)}")
```

We also provide an optimized graph data structure for partitioning plain graphs. The following example loads and partitions a graph:

```py
# Load graph from file (assumes Metis file format per default)
graph = mtkahypar.graph_from_file("path/to/graph/file", context)

# Partition graph
partitioned_graph = graph.partition(context)
```

**Note** that we internally use different data structures to represent a (hyper)graph based on the corresponding configuration (`PresetType`).
Therefore, you can not partition a (hyper)graph with all available configurations once it is loaded or constructed. However, you can check the compatibility of a hypergraph with a configuration with the following code:

```py
context = mtk.context_from_preset(mtkahypar.PresetType.QUALITY)
# Check if the hypergraph is compatible with the QUALITY preset
if hypergraph.is_compatible(context.preset):
   partitioned_hg = hypergraph.partition(context)
```

Supported Objective Functions
-----------

Mt-KaHyPar can optimize several objective functions which we explain in the following in more detail.

### Cut-Net Metric

![cut_net](https://github.com/kahypar/mt-kahypar/assets/9654047/bc7fc7c7-8ac4-4711-8aec-d0526ef2452c)

The cut-net metric is defined as total weight of all nets spanning more than one block of the partition Π (also called *cut nets*).


### Connectivity Metric

![connectivity](https://github.com/kahypar/mt-kahypar/assets/9654047/1c586ff4-63c3-4260-9ef5-98a76578be46)

The connectivity metric additionally multiplies the weight of each cut net with the number of blocks λ(e) spanned by that net minus one.
Thus, the connectivity metric tries to minimize the number of blocks connected by each net.


### Sum-of-external-Degrees Metric

![soed](https://github.com/kahypar/mt-kahypar/assets/9654047/4006fb4c-ac85-452e-a0d9-93d4dc7842ad)

The sum-of-external-degrees metric is similar to the connectivity metric, but does not subtract one from the number of blocks λ(e) spanned by a net.
A peculiarity of this objective function is that removing a net from the cut reduces the metric by 2ω(e), while reducing the connectivity by one reduces the metric only by ω(e).
Thus, the objective function prefers removing nets from the cut, while as a secondary criterion, it tries to reduce the connectivity of the nets.

### Steiner Tree Metric

![steiner_tree](https://github.com/kahypar/mt-kahypar/assets/9654047/926ef7d7-bb6b-4959-af0c-75ebd6f6299f)

The Steiner tree metric is the most versatile metric that we provide at the moment. A Steiner tree is a tree with minimal weight that connects a subset of the nodes on a graph (a more detailed definition can be found [here][SteinerTrees]).
For a subset with exactly two nodes, finding a Steiner tree reverts to computing the shortest path between the two nodes.
When optimizing the Steiner tree metric, we map the node set of a hypergraph H onto the nodes of a target graph G.
The objective is to minimize the total weight of all Steiner trees induced by the nets of H on G.
For a net e, dist(Λ(e)) is the weight of the minimal Steiner tree connecting the blocks Λ(e) spanned by net e on G.
The Steiner tree metric can be used to accurately model wire-lengths in VLSI design or communication costs in distributed systems when some processors do not communicate with each other directly or with different speeds.

Note that finding a Steiner tree is an NP-hard problem. We therefore enforce a strict upper bound on the number of nodes of the target graph G, which is 64 nodes at the moment.
If you want to map a hypergraph onto larger target graphs, you can use recursive partitioning.
For example, if you want to map a hypergraph onto a graph with 4096 nodes, you can first partition the hypergraph into 64 blocks, and then map each block of the partition onto a subgraph of the target graph with 64 nodes.

### Custom Objective Functions

Our implementation uses a common interface for all gain computation techniques that we use in our refinement algorithms.
This enables us to extend Mt-KaHyPar with new objective functions without having to modify the internal implementation of the refinement algorithms.
A step-by-step guide on how you can implement your own objective function can be found [here][CustomObjectiveFunction].

Improving Compile Times
-----------

Mt-KaHyPar implements several graph and hypergraph data structures, and supports different objective functions.
In the hot parts of the algorithm, each combination of (hyper)graph data structure and objective function is compiled separately, which notably increases the compile time.
We therefore provide the cmake preset `minimal` to disable some of the features of Mt-KaHyPar for faster compilation.

With this, only the `deterministic`, `default`, and `quality` configurations are available in combination with the cut-net or connectivity metric.
Using a disabled feature will throw an error. Note that you can only disable the features in our binary, not in the C and Python interface.

For more fine-grained control, you can directly use the corresponding cmake flags:
```cmake
-DKAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES=On/Off # enables/disables graph partitioning features
-DKAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES=On/Off # enables/disables our highest-quality configuration
-DKAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES=On/Off # enables/distables large k partitioning features
-DKAHYPAR_ENABLE_SOED_METRIC=On/Off # enables/disables sum-of-external-degrees metric
-DKAHYPAR_ENABLE_STEINER_TREE_METRIC=On/Off # enables/disables Steiner tree metric
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
[CustomObjectiveFunction]: https://github.com/kahypar/mt-kahypar/tree/master/mt-kahypar/partition/refinement/gains
[SteinerTrees]: https://en.wikipedia.org/wiki/Steiner_tree_problem#Steiner_tree_in_graphs_and_variants

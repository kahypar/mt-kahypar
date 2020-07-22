<h1 align="center">MT-KaHyPar - Multi-Threaded Karlsruhe Hypergraph Partitioning</h1>

License|Linux Build|Code Coverage
:--:|:--:|:--:
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)|[![Travis-CI Status](https://travis-ci.com/kittobi1992/mt-kahypar.svg?token=cKsYKySTDzC4fU7qcsEK&branch=master)](https://travis-ci.com/kittobi1992/mt-kahypar)|[![codecov](https://codecov.io/gh/kittobi1992/mt-kahypar/branch/master/graph/badge.svg?token=sNWRRtXZjI)](https://codecov.io/gh/kittobi1992/mt-kahypar)

What is a Hypergraph? What is Hypergraph Partitioning?
-----------
Hypergraphs are a generalization of graphs, where each (hyper)edge (also called net) can
connect more than two vertices. The *k*-way hypergraph partitioning problem is the generalization of the well-known graph partitioning problem: partition the vertex set into *k* disjoint
blocks of bounded size (at most 1 + ε times the average block size), while minimizing an
objective function defined on the nets.

The two most prominent objective functions are the cut-net and the connectivity (or λ − 1)
metrics. Cut-net is a straightforward generalization of the edge-cut objective in graph partitioning
(i.e., minimizing the sum of the weights of those nets that connect more than one block). The
connectivity metric additionally takes into account the actual number λ of blocks connected by a
net. By summing the (λ − 1)-values of all nets, one accurately models the total communication
volume of parallel sparse matrix-vector multiplication and once more gets a metric that reverts
to edge-cut for plain graphs.


<img src="https://cloud.githubusercontent.com/assets/484403/25314222/3a3bdbda-2840-11e7-9961-3bbc59b59177.png" alt="alt text" width="50%" height="50%"><img src="https://cloud.githubusercontent.com/assets/484403/25314225/3e061e42-2840-11e7-860c-028a345d1641.png" alt="alt text" width="50%" height="50%">

Requirements
-----------

The Multi-Threaded Karlsruhe Hypergraph Partitioning Framework requires:

  - A 64-bit Linux operating system.
  - A modern, ![C++14](https://img.shields.io/badge/C++-17-blue.svg?style=flat)-ready compiler such as `g++` version 9 or higher or `clang` version 11.0.3 or higher.
 - The [cmake][cmake] build system (>= 3.10).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).
 - The [Intel Thread Building Blocks][tbb] library (TBB)
 - The [Portable Hardware Locality][hwloc] library (hwloc)

The following command will install most of the required dependencies:

    sudo apt-get install libboost-program-options-dev libnuma-dev numactl libhwloc-dev moreutils linux-tools-common linux-tools-generic libtbb-dev

Building MT-KaHyPar
-----------

1. Clone the repository including submodules:

   ```git clone --depth=1 --recursive git@github.com:kittobi1992/mt-kahypar.git```
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE`
4. Run make: `make KaHyPar`

The binary will be located in `build/mt-kahypar/application/`. Per default, MT-KaHyPar uses 32-bit vertex and hyperedge IDs. If you want to partition hypergraphs with more than 4.294.967.295 vertices or hyperedges, add option `-DKAHYPAR_USE_64_BIT_IDS=ON` to the `cmake` build command.

Running MT-KaHyPar
-----------

MT-KaHyPar has several configuration parameters. We recommend to use one of the two presets which are placed in folder `config`. We provide a quality (`quality_preset.ini`) and speed preset (`speed_preset.ini`). If you want to change parameters manually, please run `./KaHyPar --help` for a detailed description of the different program options. We use the [hMetis format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf) for the input hypergraph file as well as the partition output file.

To run Mt-KaHyPar, you can use the following command:

    ./KaHyPar -h <path-to-hgr> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1 -m direct -p <path-to-config>

The partition output file will be placed in the same folder than the input hypergraph file. If you want to change the default partition output folder, add command line parameter `--partition-output-folder=path/to/folder`. Further, there are several useful options that can provide you with additional insight during and after the partitioning process:
- `--show-detailed-timings=true`: Shows detailed subtimings of each multilevel phase at the end of the partitioning process
- `--show-memory-consumption=true`: Gives detailed information on how much memory was allocated and how memory is reused throughout the algorithm
- `--enable-progress-bar=true`: Shows a progess bar during the coarsening and refinement phase

[cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
[tbb]: https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html
[hwloc]: https://www.open-mpi.org/projects/hwloc/
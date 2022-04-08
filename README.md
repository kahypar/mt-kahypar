<h1 align="center">Mt-KaHyPar - Multi-Threaded Karlsruhe Hypergraph Partitioner</h1>

License|Linux Build|Code Coverage|Code Quality
:--:|:--:|:--:|:--:
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)|[![Build Status](https://github.com/kahypar/mt-kahypar/actions/workflows/mt_kahypar_ci.yml/badge.svg)](https://github.com/kahypar/mt-kahypar/actions/workflows/mt_kahypar_ci.yml)|[![codecov](https://codecov.io/gh/kahypar/mt-kahypar/branch/master/graph/badge.svg?token=sNWRRtXZjI)](https://codecov.io/gh/kahypar/mt-kahypar)|[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/kittobi1992/mt-kahypar.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kittobi1992/mt-kahypar/context:cpp)

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

What is Mt-KaHyPar?
-----------
Mt-KaHyPar is a shared-memory multilevel hypergraph partitioning framework
for optimizing the (λ − 1)-metric.
As a multilevel algorithm, it consist of three phases: In the *coarsening phase*, the
hypergraph is coarsened to obtain a hierarchy of smaller hypergraphs. After applying an
*initial partitioning* algorithm to the smallest hypergraph in the second phase, coarsening is
undone and, at each level, several *local search* methods are used to improve the partition induced by
the coarser level. Additionally, we use a hypergraph clustering algorithm as preprocessing
to restrict contractions to densely coupled regions during coarsening.

The Mt-KaHyPar framework provides two hypergraph partitioners and a graph partitioner:

- **Mt-KaHyPar-D**: A scalable partitioner that computes good partitions very fast (for hypergraphs)
- **Mt-KaHyPar-Graph**: A scalable partitioner that computes good partitions very fast (for graphs)
- **Mt-KaHyPar-Q**: A scalable partitioner that computes high-quality partitions

Requirements
-----------

The Multi-Threaded Karlsruhe Hypergraph Partitioning Framework requires:

  - A 64-bit Linux operating system.
  - A modern, ![C++17](https://img.shields.io/badge/C++-17-blue.svg?style=flat)-ready compiler such as `g++` version 7 or higher or `clang` version 11.0.3 or higher.
 - The [cmake][cmake] build system (>= 3.16).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).
 - The [Intel Thread Building Blocks][tbb] library (TBB)
 - The [Portable Hardware Locality][hwloc] library (hwloc)

The following command will install most of the required dependencies on a Ubuntu machine:

    sudo apt-get install libboost-program-options-dev libnuma-dev numactl libhwloc-dev moreutils linux-tools-common linux-tools-generic libtbb-dev

Building Mt-KaHyPar
-----------

1. Clone the repository including submodules:

   ```git clone --depth=1 --recursive git@github.com:kahypar/mt-kahypar.git```
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE`
4. Run make: `make MtKaHyPar -j`

The build produces five executables, which will be located in `build/mt-kahypar/application/`:

- `MtKaHyParDefault` and `MtKaHyParGraph` (Mt-KaHyPar-D): computes good partitions very fast
- `MtKaHyPar(Graph)Quality` (Mt-KaHyPar-Q): computes high-quality partitions in reasonable time (using n levels)
- `MtKaHyPar`: wrapper around the five binaries

Note that `MtKaHyParGraph` and `MtKaHyParGraphQuality` uses the same feature set as `MtKaHyParDefault` and `MtKaHyParQuality`. However, they replace the internal hypergraph data structure of with a graph data structure. In fact, both are a factor of 2 faster for graphs on average.

Running Mt-KaHyPar
-----------

Mt-KaHyPar has several configuration parameters. We recommend to use one of our presets (also located in the `config` folder):

- `default`: default parameters for Mt-KaHyPar-D/-Graph (`config/default_preset.ini`)
- `default_flows`: extends the default preset with flow-based refinement (`config/default_flow_preset.ini`)
- `deterministic`: parameters to make Mt-KaHyPar-D deterministic (`config/deterministic_preset.ini`)
- `quality`: default parameters for Mt-KaHyPar-Q (`config/quality_preset.ini`)
- `quality_flows`: extends the quality preset with flow-based refinement (`config/quality_flow_preset.ini`)

The presets can be ranked from lowest to the highest quality as follows: `deterministic`,
`default`, `quality`, `default_flows` and `quality_flows`.
Deterministic mode is only supported for Mt-KaHyPar-D, not -Graph or -Q.
If you want to change parameters manually, please run `--help` for a detailed description of the different program options. We use the [hMetis format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf) for hypergraph files as well as the partition output file and the [Metis format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf) for graph files. Per default, we expect the input to be in hMetis format, but you can read graphs in Metis format via command line parameter `--input-file-format=metis`.

To run Mt-KaHyPar, you can use the following command:

    ./MtKaHyPar -h <path-to-hgr> --preset-type=<deterministic/default/default_flows/quality/quality_flows> --instance_type=<hypergraph/graph> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1 -m direct

or directly provide a configuration file (see `config` folder):

    ./MtKaHyPar -h <path-to-hgr> -p <path-to-config-file> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1 -m direct

Note that when `--instance-type=graph` is set, we run Mt-KaHyPar-Graph (only available for preset types `default` and `default_flows`), otherwise Mt-KaHyPar-D or -Q based on the preset type. The partition output file will be placed in the same folder as the input hypergraph file. If you want to change the default partition output folder, add the command line parameter `--partition-output-folder=path/to/folder`. There is also an option to disable writing the partition file `--write-partition-file=false`. Further, there are several useful options that can provide you with additional insights during and after the partitioning process:
- `--show-detailed-timings=true`: Shows detailed subtimings of each multilevel phase at the end of the partitioning process
- `--show-memory-consumption=true`: Gives detailed information on how much memory was allocated and how memory is reused throughout the algorithm
- `--enable-progress-bar=true`: Shows a progess bar during the coarsening and refinement phase

Mt-KaHyPar uses 32-bit vertex and hyperedge IDs. If you want to partition hypergraphs with more than 4.294.967.295 vertices or hyperedges, add option `-DKAHYPAR_USE_64_BIT_IDS=ON` to the `cmake` build command.

Performance
-----------

We have summarized our experimental results on an [external webpage][ExperimentalResults]. The resource provides a detailed
overview of Mt-KaHyPar's performance compared to other prominent state-of-the-art systems in terms of running time
and quality.


Using the Library Interfaces
-----------

We provide a simple C-style interface to use Mt-KaHyPar as a library.  The library can be built and installed via

```sh
make install.mtkahypar # use sudo to install system-wide
```

Note: When installing locally, the build will exit with an error due to missing permissions.
However, the library is still built successfully and is available in `lib/`.

It can be used like this:

```cpp
#include <memory>
#include <vector>
#include <iostream>

#include <libkahypar.h>

int main(int argc, char* argv[]) {

  // Initialize thread pool with 8 threads and NUMA allocation policy INTERLEAVED
  mt_kahypar_initialize_thread_pool(8, true /* activate interleaved NUMA allocation policy */ );

  // Load context from file
  mt_kahypar_context_t* context = mt_kahypar_context_new();
  mt_kahypar_configure_context_from_file(context, "path/to/config/file");

  // Setup Hypergraph
  const mt_kahypar_hypernode_id_t num_vertices = 7;
  const mt_kahypar_hyperedge_id_t num_hyperedges = 4;

  std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> hyperedge_weights =
    std::make_unique<mt_kahypar_hyperedge_weight_t[]>(4);

  // force the cut to contain hyperedge 0 and 2
  hyperedge_weights[0] = 1;  hyperedge_weights[1] = 1000;
  hyperedge_weights[2] = 1;  hyperedge_weights[3] = 1000;

  std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);

  hyperedge_indices[0] = 0; hyperedge_indices[1] = 2;
  hyperedge_indices[2] = 6; hyperedge_indices[3] = 9;
  hyperedge_indices[4] = 12;

  std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(12);

  // hypergraph from hMetis manual page 14
  hyperedges[0] = 0;  hyperedges[1] = 2;
  hyperedges[2] = 0;  hyperedges[3] = 1;
  hyperedges[4] = 3;  hyperedges[5] = 4;
  hyperedges[6] = 3;  hyperedges[7] = 4;
  hyperedges[8] = 6;  hyperedges[9] = 2;
  hyperedges[10] = 5; hyperedges[11] = 6;

  const double imbalance = 0.03;
  const mt_kahypar_partition_id_t k = 2;

  mt_kahypar_hyperedge_weight_t objective = 0;

  std::vector<mt_kahypar_partition_id_t> partition(num_vertices, -1);

  // Partition Hypergraph
  mt_kahypar_partition(num_vertices, num_hyperedges,
       	               imbalance, k, 0 /* seed */,
               	       nullptr /* unit vertex_weights */, hyperedge_weights.get(),
               	       hyperedge_indices.get(), hyperedges.get(),
       	               &objective, context, partition.data(),
                       false /* verbose output */ );

  // Print objective and block of each vertex
  std::cout << "Objective: " << objective << std::endl;
  for ( int i = 0; i != num_vertices; ++i ) {
    std::cout << "Vertex " << i << " = " << partition[i] << std::endl;
  }

  mt_kahypar_context_free(context);
}
```

If you want to load a hypergraph from a file, you can use the following code snippet:

```cpp
mt_kahypar_hypernode_id_t num_vertices = 0;
mt_kahypar_hyperedge_id_t num_hyperedges = 0;
size_t* hyperedge_indices(nullptr);
mt_kahypar_hyperedge_id_t* hyperedges(nullptr);
mt_kahypar_hypernode_weight_t* hypernode_weights(nullptr);
mt_kahypar_hyperedge_weight_t* hyperedge_weights(nullptr);
mt_kahypar_read_hypergraph_from_file("path/to/hypergraph/file", &num_vertices, &num_hyperedges,
  &hyperedge_indices, &hyperedges, &hyperedge_weights, &hypernode_weights);
```

To compile the program using `g++` and our default hypergraph partitioner (Mt-KaHyPar-D) run:

```sh
g++ -std=c++17 -DNDEBUG -O3 your_program.cc -o your_program -lmtkahypard
```

To compile the program using `g++` and our fast graph partitioner (Mt-KaHyPar-Graph) run:

```sh
g++ -std=c++17 -DNDEBUG -O3 your_program.cc -o your_program -lmtkahypargraph
```

To compile the program using `g++` and our strong hypergraph partitioner (Mt-KaHyPar-Q) run:

```sh
g++ -std=c++17 -DNDEBUG -O3 your_program.cc -o your_program -lmtkahyparq
```

To execute the binary, you need to ensure that the installation directory (probably `/usr/local/lib` for system-wide installation)
is included in the dynamic library path.
The path can be updated with:

```sh
LD_LIBRARY_PATH="$LD_LIBRARY_PATH;/usr/local/lib"
export LD_LIBRARY_PATH
```

To remove the library from your system use the provided uninstall target:

```sh
make uninstall-mtkahypar
```

You can find detailed examples how to use our library interface in the folder `examples`.

Bug Reports
-----------

We encourage you to report any problems with Mt-KaHyPar via the [github issue tracking system](https://github.com/kittobi1992/mt-kahypar/issues) of the project.

Licensing
---------

Mt-KaHyPar is free software provided under the GNU General Public License (GPLv3).
For more information see the [LICENSE file][LF].
We distribute this framework freely to foster the use and development of hypergraph partitioning tools.
If you use Mt-KaHyPar in an academic setting please cite the appropriate papers.
If you are interested in a commercial license, please contact me.

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

    // Deterministic Partitioning (Technical Report)
    @techreport{MT-KAHYPAR-SDET,
      title       = {Deterministic Parallel Hypergraph Partitioning},
      author      = {Lars Gottesbüren and
                     Michael Hamann},
      institution = {Karlsruhe Institute of Technology},
      year        = {2021},
      url         = {https://arxiv.org/pdf/2112.12704.pdf}
    }

    // Flow-Based Refinement (Technical Report)
    @techreport{MT-KAHYPAR-FLOWS,
      title       = {Parallel Flow-Based Hypergraph Partitioning},
      author      = {Lars Gottesb{\"{u}}ren and
                     Tobias Heuer and
                     Peter Sanders},
      institution = {Karlsruhe Institute of Technology},
      year        = {2022},
      url         = {https://arxiv.org/pdf/2201.01556.pdf}
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

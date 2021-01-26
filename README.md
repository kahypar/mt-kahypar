<h1 align="center">Mt-KaHyPar - Multi-Threaded Karlsruhe Hypergraph Partitioner</h1>

License|Linux Build|Code Coverage|Code Quality
:--:|:--:|:--:|:--:
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)|[![Travis-CI Status](https://travis-ci.org/kahypar/mt-kahypar.svg?branch=master)](https://travis-ci.org/kahypar/mt-kahypar)|[![codecov](https://codecov.io/gh/kahypar/mt-kahypar/branch/master/graph/badge.svg?token=sNWRRtXZjI)](https://codecov.io/gh/kahypar/mt-kahypar)|[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/kittobi1992/mt-kahypar.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kittobi1992/mt-kahypar/context:cpp)

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
Mt-KaHyPar is a shared-memory direct k-way multilevel hypergraph partitioning framework
for optimizing the (λ − 1)-metric.
As a multilevel algorithm, it consist of three phases: In the *coarsening phase*, the
hypergraph is coarsened to obtain a hierarchy of smaller hypergraphs. After applying an
*initial partitioning* algorithm to the smallest hypergraph in the second phase, coarsening is
undone and, at each level, several *local search* method are used to improve the partition induced by
the coarser level. Additionally, we use a hypergraph clustering algorithm as preprocessing
to restrict contractions to densly coupled regions of the hypergraph during the coarsening phase.

Scalability of Mt-KaHyPar
-----------
To evaluate speedups of Mt-KaHyPar, we used a benchmark set consisting of 94 large hypergraphs (see [Benchmark Statistics][SetB]) derived from
three different application domains (VLSI design, sparse matrices and SAT formulas). In the plot below, we summarize the speedups of Mt-KaHyPar
with varying number of threads (p = {4,16,64}) and different number of blocks (k = {2,8,16,64}) for each algorithmic component. We represent
the speedup of each instance as a point and the cumulative harmonic mean speedup over all instances with a single-threaded running time >= x
seconds with a line.

The overall harmonic mean speedup of Mt-KaHyPar is 3.4 for p = 4, 10.8 for p = 16 and 18.4 for p = 64. If we only consider instances with a
single-threaded running time >= 100s, we achieve a harmonic mean speedup of 23.5 for p = 64. For p = 4, the speedup is at least 3 on over
90% of our instances.

<img src="https://user-images.githubusercontent.com/9654047/88180477-bb86f500-cc2d-11ea-8bbe-834bc1b70142.png" alt="alt text" width="100%" height="100%">

Quality of Mt-KaHyPar
-----------

We use the [*performance profiles*](https://link.springer.com/article/10.1007/s101070100263) to compare Mt-KaHyPar to other partitioning algorithms in terms of solution quality.
  For a set of <img src="https://user-images.githubusercontent.com/484403/80751017-55f10400-8b29-11ea-9d73-be63c0727d39.jpg"/> algorithms and a benchmark set <img src="https://user-images.githubusercontent.com/484403/80751742-a452d280-8b2a-11ea-9f8d-47cbdd9cfccf.jpg"/> containing <img src="https://user-images.githubusercontent.com/484403/80751744-a452d280-8b2a-11ea-9049-5c3fb2d3cc76.jpg"/> instances, the *performance ratio* <img src="https://user-images.githubusercontent.com/484403/80751746-a4eb6900-8b2a-11ea-8413-ff75bccc2296.jpg"/> relates the cut computed by
  partitioner *p* for instance *i* to the smallest minimum cut of *all* algorithms, i.e.,
  <p align="center">
<img src="https://user-images.githubusercontent.com/484403/80750749-f09d1300-8b28-11ea-859d-e3b72f543ed1.png"/>. </p>

The *performance profile* <img src="https://user-images.githubusercontent.com/484403/80751752-a583ff80-8b2a-11ea-8f67-b88625b9e958.jpg"/> of algorithm *p* is then given by the function
  <p align="center">
<img src="https://user-images.githubusercontent.com/484403/80750914-35c14500-8b29-11ea-8e2c-3203b1776a96.jpg"/>.</p>

For connectivity optimization, the performance ratios are computed using the connectivity values <img src="https://user-images.githubusercontent.com/484403/80751741-a3ba3c00-8b2a-11ea-9509-6aafec2ca490.jpg"/> instead of the cut values.
    The value of <img src="https://user-images.githubusercontent.com/484403/80751750-a583ff80-8b2a-11ea-8782-3d44ee478d54.png"/> corresponds to the fraction of instances for which partitioner *p* computed the best solution, while  <img src="https://user-images.githubusercontent.com/484403/80751750-a583ff80-8b2a-11ea-8782-3d44ee478d54.png"/> is the probability
    that a performance ratio <img src="https://user-images.githubusercontent.com/484403/80751746-a4eb6900-8b2a-11ea-8413-ff75bccc2296.jpg"/> is within a factor of <img src="https://user-images.githubusercontent.com/484403/80752228-70c47800-8b2b-11ea-99e5-4524298c8620.jpg"/> of the best possible ratio.
    Note that since performance profiles only allow to assess the performance of each algorithm relative to the *best* algorithm, the <img src="https://user-images.githubusercontent.com/484403/80751747-a4eb6900-8b2a-11ea-846f-1265085ad086.jpg"/> values
    cannot be used to rank algorithms (i.e., to determine which algorithm is the second best etc.).

In our experimental analysis, the performance profile plots are based on the *best* solutions (i.e., *minimum* connectivity/cut) each algorithm found for each instance.
    Furthermore, we choose parameters <img src="https://user-images.githubusercontent.com/484403/80751754-a61c9600-8b2a-11ea-8fdb-3ba461bfe626.jpg"/> for all *p*, *i*, and <img src="https://user-images.githubusercontent.com/484403/80751762-a6b52c80-8b2a-11ea-9f6c-8fb9d00130d3.jpg"/> such that a performance ratio <img src="https://user-images.githubusercontent.com/484403/80751756-a61c9600-8b2a-11ea-88fe-22e2bdd511aa.jpg"/> if and only if algorithm *p* computed an infeasible solution
    for instance *i*, and <img src="https://user-images.githubusercontent.com/484403/80751759-a6b52c80-8b2a-11ea-9cbc-527965565037.jpg"/> if and only if the algorithm could not compute a solution for instance *i* within the given time limit. In our performance profile plots, performance ratios corresponding to *infeasible* solutions will be shown on the x-tick with a x-symbol, while
instances that could not be partitioned within the time limit are shown with a clock symbol.

To compare us against different sequential hypergraph partitioner, we use a benchmark set consisting of 488 hypergraphs (see [Benchmark Statistics][SetA]). In the figures, we compare Mt-KaHyPar with the sequential hypergraph partitioners
PaToH 3.3 in quality (PaToH-Q), default (PaToH-D) and speed mode (PaToH-S), the recursive bisection variant (hMETIS-R) of hMETIS 2.0 and
KaHyPar-CA (similiar algorithmic components as Mt-KaHyPar) and KaHyPar-HFC (extends KaHyPar-CA with flow-based refinement) of the
[KaHyPar](https://kahypar.org/) framework. On the same benchmark set on which we performed our scalability experiments
with 94 large hypergraph (see [Benchmark Statistics][SetB]), we compare us against the distributed hypergraph partitioner Zoltan 3.83 and the default mode of PaToH 3.3.

Note that increasing number of threads does not negatively affect solution quality of Mt-KaHyPar.

<img src="https://user-images.githubusercontent.com/9654047/88182577-a6f82c00-cc30-11ea-9234-752149cab2a9.png" alt="alt text" width="100%" height="100%">
<img src="https://user-images.githubusercontent.com/9654047/88186573-a1e9ab80-cc35-11ea-8e43-e226dadc65cc.png" alt="alt text" width="100%" height="100%">

Requirements
-----------

The Multi-Threaded Karlsruhe Hypergraph Partitioning Framework requires:

  - A 64-bit Linux operating system.
  - A modern, ![C++14](https://img.shields.io/badge/C++-17-blue.svg?style=flat)-ready compiler such as `g++` version 7 or higher or `clang` version 11.0.3 or higher.
 - The [cmake][cmake] build system (>= 3.16).
 - The [Boost - Program Options][Boost.Program_options] library and the boost header files (>= 1.48).
 - The [Intel Thread Building Blocks][tbb] library (TBB)
 - The [Portable Hardware Locality][hwloc] library (hwloc)

The following command will install most of the required dependencies:

    sudo apt-get install libboost-program-options-dev libnuma-dev numactl libhwloc-dev moreutils linux-tools-common linux-tools-generic libtbb-dev

Building Mt-KaHyPar
-----------

1. Clone the repository including submodules:

   ```git clone --depth=1 --recursive git@github.com:kittobi1992/mt-kahypar.git```
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE`
4. Run make: `make MtKaHyPar -j8`

The build produces two executables, which will be located in `build/mt-kahypar/application/`:

- `MtKaHyParFast`: A scalable hypergraph partitioner that computes good partitions very fast
- `MtKaHyParStrong`: A scalable hypergraph partitioner that computes high-quality partitions

Per default, Mt-KaHyPar uses 32-bit vertex and hyperedge IDs. If you want to partition hypergraphs with more than 4.294.967.295 vertices or hyperedges, add option `-DKAHYPAR_USE_64_BIT_IDS=ON` to the `cmake` build command.

Running Mt-KaHyPar
-----------

Mt-KaHyPar has several configuration parameters. We recommend to use our presets which are placed in the `config` folder:

- `fast_preset.ini`: Contains the default parameters for Mt-KaHyPar Fast (`MtKaHyParFast`)
- `strong_preset.ini`: Contains the default parameters for Mt-KaHyPar Strong (`MtKaHyParStrong`)

If you want to change parameters manually, please run `./MtKaHyParFast --help` or `./MtKaHyParStrong --help` for a detailed description of the different program options. We use the [hMetis format](http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/manual.pdf) for the input hypergraph file as well as the partition output file.

To run Mt-KaHyPar Fast, you can use the following command:

    ./MtKaHyParFast -h <path-to-hgr> -p <path to fast_preset.ini> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1 -m direct

To run Mt-KaHyPar Strong, you can use the following command:

    ./MtKaHyParStrong -h <path-to-hgr> -p <path to strong_preset.ini> -t <# threads> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1 -m direct

The partition output file will be placed in the same folder than the input hypergraph file. If you want to change the default partition output folder, add command line parameter `--partition-output-folder=path/to/folder`. Further, there are several useful options that can provide you with additional insights during and after the partitioning process:
- `--show-detailed-timings=true`: Shows detailed subtimings of each multilevel phase at the end of the partitioning process
- `--show-memory-consumption=true`: Gives detailed information on how much memory was allocated and how memory is reused throughout the algorithm
- `--enable-progress-bar=true`: Shows a progess bar during the coarsening and refinement phase

Using the Library Interfaces
-----------

We provide a simple C-style interface to use Mt-KaHyPar as a library.  The library can be built and installed via

```sh
make install.library
```

and can be used like this:

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

To compile the program using `g++` run:

```sh
g++ -std=c++17 -DNDEBUG -O3 your_program.cc -o your_program -lkahypar
```

To remove the library from your system use the provided uninstall target:

```sh
make uninstall-kahypar
```

Note, our library interfaces uses Mt-KaHyPar Fast. We are currently working on a solution to also integrate Mt-KaHyPar Strong.

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

    // Mt-KaHyPar Fast
    @inproceedings{MT-KAHYPAR-FAST,
      title     = {Scalable Shared-Memory Hypergraph Partitioning},
      author    = {Gottesbüren, Lars and
                   Heuer, Tobias and
                   Sanders, Peter and
                   Schlag, Sebastian},
      booktitle = {23rd Workshop on Algorithm Engineering and Experiments, (ALENEX 2021)},
      pages     = {16--30},
      year      = {2021},
      publisher = {SIAM}
    }

Contributing
------------
If you are interested in contributing to the Mt-KaHyPar framework
feel free to contact me or create an issue on the
[issue tracking system](https://github.com/kittobi1992/mt-kahypar/issues).

[cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
[tbb]: https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html
[hwloc]: https://www.open-mpi.org/projects/hwloc/
[LF]: https://github.com/kittobi1992/mt-kahypar/blob/master/LICENSE "Licence"
[SetA]: http://algo2.iti.kit.edu/heuer/alenex21/instances.html?benchmark=set_a
[SetB]: http://algo2.iti.kit.edu/heuer/alenex21/instances.html?benchmark=set_b

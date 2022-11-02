/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "mt-kahypar/partition/initial_partitioning/deep_initial_partitioner.h"

#include <algorithm>
#include <limits>
#include <vector>

#include "tbb/parallel_invoke.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/coarsening/multilevel_uncoarsener.h"
#include "mt-kahypar/partition/coarsening/nlevel_uncoarsener.h"
#include "mt-kahypar/partition/deep_multilevel.h"

#include "mt-kahypar/partition/initial_partitioning/flat/pool_initial_partitioner.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"


namespace mt_kahypar {
  DeepInitialPartitioner::DeepInitialPartitioner(PartitionedHypergraph& hypergraph,
                                                 const Context& context) :
    _hg(hypergraph),
    _context(context) { }

  void DeepInitialPartitioner::initialPartitionImpl() {
    deep_multilevel::partition(_hg, _context);
  }
} // namepace mt_kahypar

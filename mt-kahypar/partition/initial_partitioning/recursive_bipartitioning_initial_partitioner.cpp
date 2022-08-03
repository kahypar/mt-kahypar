/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/initial_partitioning/recursive_bipartitioning_initial_partitioner.h"


#include <algorithm>
#include <vector>

#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/separated_nodes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/multilevel.h"

#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/recursive_bipartitioning.h"
#include "mt-kahypar/partition/star_partitioning/star_partitioning.h"

namespace mt_kahypar {
  using ds::SeparatedNodes;
  using ds::Array;

  RecursiveBipartitioningInitialPartitioner::RecursiveBipartitioningInitialPartitioner(PartitionedHypergraph& hypergraph,
                                                                             const Context& context) :
    _hg(hypergraph),
    _context(context),
    _original_s_nodes(&hypergraph.separatedNodes()),
    _s_nodes(hypergraph.separatedNodes().coarsest().copy(parallel_tag_t())) {
    }

  void RecursiveBipartitioningInitialPartitioner::initialPartitionImpl() {
    _s_nodes.onliest().cleanBatchState();
    _s_nodes.onliest().setSavepoint();
    _hg.setSeparatedNodes(&_s_nodes);
    recursive_bipartitioning::partition(_hg, _context);
    _hg.setSeparatedNodes(_original_s_nodes);
    _hg.resetSeparatedParts();

    // TODO: compare quality?!
    star_partitioning::partition(_hg, _context);
    _hg.updateBlockWeights(); // TODO: not necessary anymore?
  }
} // namepace mt_kahypar
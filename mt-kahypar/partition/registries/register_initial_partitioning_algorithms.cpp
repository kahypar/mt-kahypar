/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"

#include "mt-kahypar/partition/initial_partitioning/deep_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/recursive_bipartitioning_initial_partitioner.h"

#define REGISTER_INITIAL_PARTITIONER(id, partitioner)                                           \
  static kahypar::meta::Registrar<InitialPartitionerFactory> register_ ## partitioner(          \
    id,                                                                                         \
    [](PartitionedHypergraph& hypergraph, const Context& context)                               \
    -> IInitialPartitioner* {                                                                   \
    return new partitioner(hypergraph, context);                                                \
  })

namespace mt_kahypar {
REGISTER_INITIAL_PARTITIONER(Mode::deep_multilevel, DeepInitialPartitioner);
REGISTER_INITIAL_PARTITIONER(Mode::recursive_bipartitioning, RecursiveBipartitioningInitialPartitioner);
}  // namespace mt_kahypar

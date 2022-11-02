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
#include "mt-kahypar/partition/preprocessing/sparsification/hypergraph_sparsifier.h"
#include "mt-kahypar/partition/preprocessing/sparsification/policies/similiar_net_combine.h"

#define REGISTER_HYPERGRAPH_SPARSIFIER(id, sparsifier)                                      \
  static kahypar::meta::Registrar<HypergraphSparsifierFactory> register_ ## sparsifier(     \
    id,                                                                                     \
    [](const Context& context)                                                              \
    -> IHypergraphSparsifier* {                                                             \
    return new sparsifier(context);                                                         \
  })

namespace mt_kahypar {
using HypergraphUnionSparsifier = HypergraphSparsifier<UnionCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::union_nets, HypergraphUnionSparsifier);
using HypergraphMaxSizeSparsifier = HypergraphSparsifier<MaxSizeCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::max_size, HypergraphMaxSizeSparsifier);
using HypergraphImportanceSparsifier = HypergraphSparsifier<NetImportanceCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::importance, HypergraphImportanceSparsifier);
using HypergraphUndefinedSparsifier = HypergraphSparsifier<UndefinedCombiner>;
REGISTER_HYPERGRAPH_SPARSIFIER(SimiliarNetCombinerStrategy::UNDEFINED, HypergraphUndefinedSparsifier);
}  // namespace mt_kahypar

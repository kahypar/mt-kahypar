/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

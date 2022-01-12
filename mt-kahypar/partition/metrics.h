/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar::metrics {

HyperedgeWeight hyperedgeCut(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight km1(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight soed(const PartitionedHypergraph& hypergraph, bool parallel = true);

bool isBalanced(const PartitionedHypergraph& phg, const Context& context);

HyperedgeWeight objective(
        const PartitionedHypergraph& hg,
        const kahypar::Objective& objective,
        bool parallel = true);

double imbalance(const PartitionedHypergraph& hypergraph, const Context& context);

HyperedgeWeight judiciousLoad(const PartitionedHypergraph& hypergraph, bool parallel = true);

HyperedgeWeight minLoad(const PartitionedHypergraph& hypergraph, const bool parallel = true);

HyperedgeID maxHnDegree(const Hypergraph& hypergraph);

}  // namespace mt_kahypar

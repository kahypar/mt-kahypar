/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2022 Nikolai Maas <nikolai.maas@student.kit.edu>
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

#include "rate_separated.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/star_partitioning/star_partitioning.h"

namespace mt_kahypar {
namespace star_partitioning {

HyperedgeWeight rateSeparated(PartitionedHypergraph& phg, const Context& context, const HypernodeID u, const PartitionID to) {
    const PartitionID from = phg.partID(u);
    ASSERT(context.initial_partitioning.rater == IPSNodesRater::naive);
    if (from == to) {
        return 0;
    }
    const HyperedgeID added_old = star_partitioning::partition(phg, phg.separatedNodes().finest(), context, false);
    phg.resetSeparatedParts(false);
    if (from == kInvalidPartition) {
        phg.setNodePart(u, to);
    } else {
        phg.changeNodePart(u, from, to);
    }
    const HyperedgeID added_new = star_partitioning::partition(phg, phg.separatedNodes().finest(), context, false);
    phg.resetSeparatedParts(false);
    if (from == kInvalidPartition) {
        phg.resetPart(u, to);
    } else {
        phg.changeNodePart(u, to, from);
    }
    ASSERT(phg.partID(u) == from);
    return added_old - added_new;
}

} // namepace star_partitioning
} // namepace mt_kahypar

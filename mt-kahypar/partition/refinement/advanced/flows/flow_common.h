/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "datastructure/flow_hypergraph_builder.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {

struct FlowProblem {
  whfc::Node source;
  whfc::Node sink;
  HyperedgeWeight total_cut;
  HyperedgeWeight non_removable_cut;
  HypernodeWeight weight_of_block_0;
  HypernodeWeight weight_of_block_1;
};


} // namespace mt_kahypar
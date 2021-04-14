/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "deterministic_multilevel_coarsener.h"

#include "mt-kahypar/utils/progress_bar.h"

namespace mt_kahypar {

void DeterministicMultilevelCoarsener::coarsenImpl() {
  HypernodeID initial_num_nodes = Base::currentNumNodes();
  utils::ProgressBar progress_bar(initial_num_nodes, 0,
                                  _context.partition.verbose_output && _context.partition.enable_progress_bar);

  progress_bar += (initial_num_nodes - progress_bar.count());
  progress_bar.disable();
}

}

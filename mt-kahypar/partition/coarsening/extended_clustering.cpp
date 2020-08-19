/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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


#include "extended_clustering.h"
#include <iostream>

namespace mt_kahypar {



PartitionedHypergraph&& ExtendedClustering::uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                      std::unique_ptr<IRefiner>& fm) {
  return Base::doUncoarsen(label_propagation, fm);
}



void ExtendedClustering::coarsenImpl() {
  std::cout << "Not implemented yet. Use multilevel_coarsener instead"; std::exit(0);
}


}
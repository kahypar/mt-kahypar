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

#pragma once


#include "mt-kahypar/partition/context.h"

namespace mt_kahypar::io {
  void printStripe();
  void printBanner();
  void printCutMatrix(const PartitionedHypergraph& hypergraph);
  void printHypergraphInfo(const Hypergraph& hypergraph,
                           const std::string& name,
                           const bool show_memory_consumption);
  void printPartitioningResults(const PartitionedHypergraph& hypergraph,
                                const Context& context,
                                const std::string& description);
  void printPartitioningResults(const PartitionedHypergraph& hypergraph,
                                const Context& context,
                                const std::chrono::duration<double>& elapsed_seconds);
  void printPartWeightsAndSizes(const PartitionedHypergraph& hypergraph, const Context& context);
  void printContext(const Context& context);
  void printMemoryPoolConsumption(const Context& context);
  void printInputInformation(const Context& context, const Hypergraph& hypergraph);
  void printTopLevelPreprocessingBanner(const Context& context);
  void printCommunityInformation(const Hypergraph& hypergraph);
  void printCoarseningBanner(const Context& context);
  void printInitialPartitioningBanner(const Context& context);
  void printLocalSearchBanner(const Context& context);
  void printVCycleBanner(const Context& context, const size_t vcycle_num);
}  // namespace mt_kahypar::io

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

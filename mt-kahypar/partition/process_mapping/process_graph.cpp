/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/process_mapping/process_graph.h"

#include <cmath>
#include <limits>

#include "mt-kahypar/partition/process_mapping/all_pair_shortest_path.h"

namespace mt_kahypar {

void ProcessGraph::precomputeDistances(const size_t max_connectivity) {
  const size_t num_entries = std::pow(_k, max_connectivity);
  if ( num_entries > MEMORY_LIMIT ) {
    ERR("Too much memory requested for precomputing distances"
      << "of connectivity set in process graph.");
  }
  _distances.assign(num_entries, std::numeric_limits<HyperedgeWeight>::max() / 3);
  AllPairShortestPath::compute(_graph, _distances);

  _is_initialized = true;
}

}  // namespace kahypar

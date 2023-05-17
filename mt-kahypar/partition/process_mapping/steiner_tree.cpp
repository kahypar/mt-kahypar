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

#include "mt-kahypar/partition/process_mapping/steiner_tree.h"
#include "mt-kahypar/partition/process_mapping/all_pair_shortest_path.h"

namespace mt_kahypar {

void SteinerTree::compute(const ds::StaticGraph& graph,
                          const size_t max_set_size,
                          vec<HyperedgeWeight>& distances) {
  unused(max_set_size);
  AllPairShortestPath::compute(graph, distances);

  /**
   * Algorithm idea:
   *
   * S = distances
   * for m = 2 to max_set_size - 1 do
   *   for all subsets D of V with |D| = m do
   *     for each u in V do
   *       min_dist = inf
   *       for each subset proper subset E of D do
   *         F = D \ E
   *         min_dist = min( min_dist, S[ E u { u } ] + S[ F u { u } ] )
   *       for each v in V do
   *         S[ D u { v } ] = min( S[ D u { v } ], S[ {u, v} ] + min_dist )
   *
   *
   */
}

}  // namespace kahypar

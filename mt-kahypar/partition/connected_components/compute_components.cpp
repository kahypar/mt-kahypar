/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
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

#include "mt-kahypar/partition/connected_components/compute_components.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace connected_components {

using Bitset = mt_kahypar::ds::Bitset;

template<typename PartitionedHypergraph>
void compute_components_per_block(const PartitionedHypergraph& phg,
                                  const Context& context,
                                  vec<vec<ConnectedComponent>>& result) {
  result.resize(phg.k());
  for (vec<ConnectedComponent>& list: result) {
    list.clear();
  }

  (void)context;

  Bitset node_colored;
  node_colored.resize(phg.initialNumNodes());
  
  std::queue<HypernodeID> node_queue;
  PartitionID current_partition;

  for (const HypernodeID& hn : phg.nodes()) {
    if (node_colored.isSet((size_t) hn)) {
      continue;
    }

    current_partition = phg.partID(hn);
    node_queue.push(hn);
    
    ConnectedComponent cc = { };

    while (node_queue.size() > 0) {
      HypernodeID current = node_queue.front();
      node_colored.set((size_t) current);
      cc.nodes.push_back(current);
      node_queue.pop();

      for (const HyperedgeID& he : phg.incidentEdges(current)) {
        for (const HypernodeID& incident_hn : phg.pins(he)) {
          if (node_colored.isSet((size_t) incident_hn)) {
            continue;
          }
          
          if (phg.partID(incident_hn) != current_partition) {
            continue;
          }

          node_queue.push(incident_hn);
        }
      }
    }
    
    result[current_partition].push_back(cc);
  }
}



namespace {
#define COMPUTE_COMPONENTS_PER_BLOCK(X) void compute_components_per_block(const X& phg, const Context& context, vec<vec<ConnectedComponent>>& result)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(COMPUTE_COMPONENTS_PER_BLOCK)

}  // namespace connected_components
}  // namespace mt_kahypar

#include "mt-kahypar/partition/k_colouring.h"
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

template<typename TypeTraits>
void KColouring<TypeTraits>::partition(PartitionedHypergraph& phg,
                                        const Context& context) {
  Hypergraph& hg = phg.hypergraph();

  const ds::DynamicGraph& constraint_graph = phg.fixedVertexSupport().getConstraintGraph();
  NodeSelector selector(constraint_graph);
  graph_colouring colouring;
  colouring.node_colours = vec<Colour> (constraint_graph.numNodes(), kInvalidColour);
  colouring.used_colours = phg.k();
  block_weights weights;
  weights.weights = vec<HypernodeWeight> (phg.k(), 0);

  for (HypernodeID hn = selector.getNextNode(); hn != kInvalidHypernode; hn = selector.getNextNode()) {
    HypernodeID hn_id = constraint_graph.nodeWeight(hn);
    vec<bool> is_colour_usable = getPossibleColours(constraint_graph.incidentNodes(hn), colouring);
    //select best partition
    Colour balanced_colour = getMostBalancedColor(is_colour_usable, weights);
    if (balanced_colour >= phg.k()) {
      LOG << "trying to put node in partition"<<balanced_colour;
      balanced_colour = balanced_colour % phg.k();
    }
    weights.setNodePart(phg, hn_id, balanced_colour);
    colouring.node_colours[hn] = balanced_colour;
    phg.setOnlyNodePart(hn_id, balanced_colour);
  }
  phg.initializePartition();
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(KColouring)

}// namespace mt_kahypar
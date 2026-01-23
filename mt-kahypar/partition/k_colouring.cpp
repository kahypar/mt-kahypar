#include "mt-kahypar/partition/k_colouring.h"
#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

template<typename TypeTraits>
void KColouring<TypeTraits>::partition() {
  const ds::DynamicGraph& constraint_graph = _phg.fixedVertexSupport().getConstraintGraph();
  NodeSelector selector(constraint_graph);
  graph_colouring colouring;
  colouring.node_colours = vec<Colour> (constraint_graph.numNodes(), kInvalidColour);
  colouring.used_colours = _phg.k();
  block_weights weights;
  weights.weights = vec<HypernodeWeight> (_phg.k(), 0);
  LOG << "max weight:"<< _context.partition.max_part_weights[0];

  for (HypernodeID hn = selector.getNextNode(); hn != kInvalidHypernode; hn = selector.getNextNode()) {
    HypernodeID hn_id = constraint_graph.nodeWeight(hn);
    vec<bool> is_colour_usable = getPossibleColours(constraint_graph.incidentNodes(hn), colouring);
    //select best partition
    Colour max_balanced_gain_colour = getBestCutColor(is_colour_usable, hn_id);
    if (max_balanced_gain_colour >= _phg.k()) {
      LOG << "trying to put node in partition"<<max_balanced_gain_colour;
      max_balanced_gain_colour = max_balanced_gain_colour % _phg.k();
    }
    weights.setNodePart(_phg, hn_id, max_balanced_gain_colour);
    colouring.node_colours[hn] = max_balanced_gain_colour;
    _phg.setNodePart(hn_id, max_balanced_gain_colour);
  }
  setRemainingNodes();
}

template<typename TypeTraits>
void KColouring<TypeTraits>::setRemainingNodes() {
  vec<bool> allColoursUsable(_phg.k(), true);
  for (const HypernodeID hn_id : _phg.nodes()) {
    if (_phg.partID(hn_id) == kInvalidPartition) {
      PartitionID bestPartition = getBestCutColor(allColoursUsable, hn_id);
      _phg.setNodePart(hn_id, bestPartition);
    }
  }
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(KColouring)

}// namespace mt_kahypar
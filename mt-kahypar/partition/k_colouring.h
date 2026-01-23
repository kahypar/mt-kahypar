#pragma once

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/dynamic_graph.h"
#include "mt-kahypar/datastructures/priority_queue.h"

namespace mt_kahypar {

using DynamicGraph = typename ds::DynamicGraph;
using Colour = int32_t;

static constexpr Colour kInvalidColour = -1;

struct graph_colouring {
  vec<Colour> node_colours;
  Colour used_colours = 0;
} ;

struct block_weights {
  vec<HypernodeWeight> weights;

  template<typename PartitionedHypergraph>
  void setNodePart(const PartitionedHypergraph& phg, const HypernodeID hn, const PartitionID part) {
    HypernodeWeight weight = phg.nodeWeight(hn);
    weights[part] += weight;
  }
};

class NodeSelector {
 public:
  
  struct Key {
    HypernodeID degree;
    HypernodeID id;

    bool operator>(const Key& other) const {
      if (degree != other.degree) {
        return degree > other.degree;
      }
      return id > other.id;
    }
    bool operator<(const Key& other) const {
      if (degree != other.degree) {
        return degree < other.degree;
      }
      return id < other.id;
    }
    bool operator==(const Key& other) const {
        return degree == other.degree && id == other.id;
    }
  };


  NodeSelector(const DynamicGraph& constraint_graph) :
    _constraint_graph(constraint_graph.copy()),
    _positions(constraint_graph.numNodes()),
    _pq(_positions.data(), constraint_graph.numNodes())
    {
      _constraint_graph.removeSinglePinAndParallelHyperedges();
      for (const HypernodeID& node : _constraint_graph.nodes()) {
        _pq.insert(node, Key{_constraint_graph.nodeDegree(node), node});
        if (_constraint_graph.nodeDegree(node) > 19) {
          LOG << "node"<<node<<"has a degree of"<<_constraint_graph.nodeDegree(node);
        }
      }
      LOG <<"Nodes in pq"<< _pq.size();
    }

  HypernodeID getNextNode() {
    if (_pq.empty()) {
      return kInvalidHypernode;
    }
    const HypernodeID node = _pq.top();
    _pq.deleteTop();
    // for (const auto& neighbour : _constraint_graph.incidentNodes(node)) {
    //   if (_pq.contains(neighbour)) {
    //     Key key = _pq.getKey(neighbour);
    //     if (key.degree > 0) {
    //       key.degree--;
    //     }
    //     _pq.adjustKey(neighbour, key);
    //   }
    // }
    return node;
  }

 private:
  DynamicGraph _constraint_graph;

  vec<PosT> _positions;
  ds::Heap<Key, HypernodeID, std::less<Key>, 4> _pq;
};

template<typename TypeTraits>
class KColouring{

  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  KColouring(const Context& context) :
    _context(context) { }

  Colour colour(const PartitionedHypergraph& phg) {
    const ds::DynamicGraph& constraint_graph = phg.fixedVertexSupport().getConstraintGraph();
    NodeSelector selector(constraint_graph);
    graph_colouring colouring;
    colouring.node_colours = vec<Colour> (constraint_graph.numNodes(), kInvalidColour);
    for (HypernodeID hn = selector.getNextNode(); hn != kInvalidHypernode; hn = selector.getNextNode()) {
      const Colour colour = getSmalestPossibleColour(constraint_graph.incidentNodes(hn), colouring);
      colouring.node_colours[hn] = colour;
      if (colour >= colouring.used_colours) {
        colouring.used_colours++;
      }
    }
    LOG <<"used Colours"<<colouring.used_colours;
    vec<HypernodeID> count(colouring.used_colours, 0);
    for (auto node_colour : colouring.node_colours) {
      count[node_colour]++;
    }
    for(int i = 0; i<colouring.used_colours; i++) {
      LOG<<"Colour"<<i<<"num Nodes"<<count[i];
    }
    return colouring.used_colours;
  }

  static void partition(PartitionedHypergraph& hypergraph,
                        const Context& context);

 private:

  static vec<bool> getPossibleColours(IteratorRange<ds::IncidentNodeIterator> incident_iterator, graph_colouring& colouring) {
    vec<bool> is_colour_usable(colouring.used_colours, true);
    for (auto incident_hn : incident_iterator) {
      Colour incident_colour = colouring.node_colours[incident_hn];
      if (incident_colour != kInvalidColour) {
        is_colour_usable[incident_colour] = false;
      }
    }
    return is_colour_usable;
  }

  static Colour getSmalestPossibleColour(IteratorRange<ds::IncidentNodeIterator> incident_iterator, graph_colouring& colouring) {
    vec<bool> is_colour_usable = getPossibleColours(incident_iterator, colouring);
    Colour colour_tu_use = 0;
    while(!is_colour_usable[colour_tu_use]) colour_tu_use++;
    return colour_tu_use;
  }

  static Colour getMostBalancedColor(const vec<bool>& is_colour_usable, const block_weights& weights) {
    Colour colour_tu_use = 0;
    HypernodeWeight best_weight = std::numeric_limits<HypernodeWeight>::max();
    for (Colour colour = 0; colour < is_colour_usable.size(); colour++) {
      if (is_colour_usable[colour] && weights.weights[colour] < best_weight) { // if no color is valid all put in 0
        colour_tu_use = colour;
        best_weight = weights.weights[colour];
      }
    }
    return colour_tu_use;
  }

  static Colour getBestCutColor(const vec<bool>& is_colour_usable, const block_weights& weights) {
    Colour colour_tu_use = 0;
    HypernodeWeight best_weight = std::numeric_limits<HypernodeWeight>::max();
    for (Colour colour = 0; colour < is_colour_usable.size(); colour++) {
      if (is_colour_usable[colour] && weights.weights[colour] < best_weight) { // if no color is valid all put in 0
        colour_tu_use = colour;
        best_weight = weights.weights[colour];
      }
    }
    return colour_tu_use;
  }

  const Context& _context;
};

}// namespace mt_kahypar
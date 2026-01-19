#pragma once

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
      }
      LOG << "PQ size:"<<_pq.size();
      LOG << "graph nodes:"<<_constraint_graph.numNodes();
    }

  HypernodeID getNextNode() {
    if (_pq.empty()) {
      return kInvalidHypernode;
    }
    const HypernodeID node = _pq.top();
    _pq.deleteTop();
    for (const auto& neighbour : _constraint_graph.incidentNodes(node)) {
      if (_pq.contains(neighbour)) {
        Key key = _pq.getKey(neighbour);
        if (key.degree > 0) {
          key.degree--;
        }
        _pq.adjustKey(neighbour, key);
      }
    }
    return node;
  }

 private:
  DynamicGraph _constraint_graph;

  vec<PosT> _positions;
  ds::Heap<Key, HypernodeID, std::greater<Key>, 4> _pq;
};

template<typename TypeTraits>
class KColouring{

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

 private:

  Colour getSmalestPossibleColour(IteratorRange<ds::IncidentNodeIterator> incident_iterator, graph_colouring& colouring) {
    vec<Colour> is_colour_used(colouring.used_colours + 1, false);
    Colour colour_tu_use = 0;
    for (auto incident_hn : incident_iterator) {
      Colour incident_colour = colouring.node_colours[incident_hn];
      if (incident_colour != kInvalidColour) {
        is_colour_used[incident_colour] = true;
        if ( incident_colour == colour_tu_use) {
          while (is_colour_used[colour_tu_use])
          {
            colour_tu_use++;
          }
        }
      }
    }
    return colour_tu_use;
  }

  const Context& _context;
};

}// namespace mt_kahypar
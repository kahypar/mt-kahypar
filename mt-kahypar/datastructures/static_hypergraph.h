/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
#include "tbb/parallel_scan.h"

#include "kahypar/meta/mandatory.h"
#include "kahypar/utils/math.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace ds {

template <typename TBBNumaArena = Mandatory>
class StaticHypergraph {

  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex.
   */
  class Hypernode {
   public:
    using IDType = HyperedgeID;

    Hypernode() :
      _begin(0),
      _size(0),
      _weight(1),
      _community_id(0),
      _valid(false) { }

    Hypernode(const size_t begin,
              const size_t size,
              const HypernodeWeight weight) :
      _begin(begin),
      _size(size),
      _weight(weight),
      _community_id(0),
      _valid(true) { }

    bool isDisabled() const {
      return _valid == false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    // ! Returns the index of the first element in _incident_nets
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incident_nets to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incident_nets
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    HyperedgeWeight weight() const {
      ASSERT(!isDisabled());
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    PartitionID communityID() const {
      ASSERT(!isDisabled());
      return _community_id;
    }

    void setCommunityID(const PartitionID community_id) {
      ASSERT(!isDisabled());
      _community_id = community_id;
    }

   private:
    // ! Index of the first element in _incident_nets
    size_t _begin;
    // ! Number of incident nets
    size_t _size;
    // ! Hypernode weight
    HyperedgeWeight _weight;
    // ! Community id
    PartitionID _community_id;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /**
   * Represents a hyperedge of the hypergraph and contains all information
   * associated with a net (except connectivity information).
   */
  class Hyperedge {
   public:

    Hyperedge() :
      _begin(0),
      _size(0),
      _weight(1),
      _hash(kEdgeHashSeed),
      _valid(false) { }

    Hyperedge(const size_t begin,
              const size_t size,
              const HyperedgeWeight weight) :
      _begin(begin),
      _size(size),
      _weight(weight),
      _hash(kEdgeHashSeed),
      _valid(true) { }

    // ! Disables the hypernode/hyperedge. Disable hypernodes/hyperedges will be skipped
    // ! when iterating over the set of all nodes/edges.
    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    bool isDisabled() const {
      return _valid == false;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incidence_array to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    HyperedgeWeight weight() const {
      ASSERT(!isDisabled());
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    size_t& hash() {
      return _hash;
    }

    size_t hash() const {
      return _hash;
    }

    bool operator== (const Hyperedge& rhs) const {
      return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
    }

    bool operator!= (const Hyperedge& rhs) const {
      return !operator== (this, rhs);
    }

   private:
    // ! Index of the first element in _incidence_array
    size_t _begin;
    // ! Number of pins
    size_t _size;
    // ! hyperedge weight
    HyperedgeWeight _weight;
    // ! Hash of pins
    size_t _hash;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  static_assert(std::is_trivially_copyable<Hypernode>::value, "Hypernode is not trivially copyable");
  static_assert(std::is_trivially_copyable<Hyperedge>::value, "Hyperedge is not trivially copyable");

  using IncidenceArray = parallel::scalable_vector<HypernodeID>;
  using IncidentNets = parallel::scalable_vector<HyperedgeID>;
  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using Counter = parallel::scalable_vector<size_t>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<Counter>;

 public:
  static constexpr bool is_static_hypergraph = true;
  static constexpr bool is_partitioned = false;

  explicit StaticHypergraph() :
    _num_hypernodes(0),
    _num_hyperedges(0),
    _num_pins(0),
    _total_weight(0),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array() { }

  explicit StaticHypergraph(const HypernodeID num_hypernodes,
                            const HyperedgeID num_hyperedges,
                            const HyperedgeVector& edge_vector,
                            const HyperedgeWeight* hyperedge_weight = nullptr,
                            const HypernodeWeight* hypernode_weight = nullptr) :
    _num_hypernodes(num_hypernodes),
    _num_hyperedges(num_hyperedges),
    _num_pins(0),
    _total_weight(0),
    _hypernodes(num_hypernodes),
    _incident_nets(),
    _hyperedges(num_hyperedges),
    _incidence_array() {
    construct(edge_vector, hyperedge_weight, hypernode_weight);
  }

  StaticHypergraph(const StaticHypergraph&) = delete;
  StaticHypergraph & operator= (const StaticHypergraph &) = delete;

  StaticHypergraph(StaticHypergraph&& other) :
    _num_hypernodes(other._num_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_pins(other._num_pins),
    _total_weight(other._total_weight),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)) { }

  StaticHypergraph & operator= (StaticHypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_pins = other._num_pins;
    _total_weight = other._total_weight;
    _hypernodes = std::move(other._hypernodes);
    _incident_nets = std::move(other._incident_nets);
    _hyperedges = std::move(other._hyperedges);
    _incidence_array = std::move(other._incidence_array);
  }

 private:

  // ####################### Construction #######################

  void construct(const HyperedgeVector& edge_vector,
                 const HyperedgeWeight* hyperedge_weight = nullptr,
                 const HypernodeWeight* hypernode_weight = nullptr) {
    ASSERT(edge_vector.size() == _num_hyperedges);

    // Compute number of pins per hyperedge and number
    // of incident nets per vertex
    Counter num_pins_per_hyperedge(_num_hyperedges, 0);
    ThreadLocalCounter local_incident_nets_per_vertex(_num_hypernodes, 0);
    tbb::parallel_for(0UL, _num_hyperedges, [&](const size_t pos) {
      Counter& num_incident_nets_per_vertex = local_incident_nets_per_vertex.local();
      num_pins_per_hyperedge[pos] = edge_vector[pos].size();
      for ( const HypernodeID& pin : edge_vector[pos] ) {
        ASSERT(pin < _num_hypernodes);
        ++num_incident_nets_per_vertex[pin];
      }
    });

    // We sum up the number of incident nets per vertex only thread local.
    // To obtain the global number of incident nets per vertex, we iterate
    // over each thread local counter and sum it up.
    Counter num_incident_nets_per_vertex(_num_hypernodes, 0);
    for ( Counter& c : local_incident_nets_per_vertex ) {
      tbb::parallel_for(0UL, _num_hypernodes, [&](const size_t pos) {
        num_incident_nets_per_vertex[pos] += c[pos];
      });
    }

    // Compute prefix sum over the number of pins per hyperedge and the
    // number of incident nets per vertex. The prefix sum is used than as
    // start position for each hyperedge resp. hypernode in the incidence
    // resp. incident nets array.
    parallel::TBBPrefixSum<size_t> pin_prefix_sum(num_pins_per_hyperedge);
    parallel::TBBPrefixSum<size_t> incident_net_prefix_sum(num_incident_nets_per_vertex);
    tbb::parallel_invoke([&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, _num_hyperedges), pin_prefix_sum);
    }, [&] {
      tbb::parallel_scan(tbb::blocked_range<size_t>(
        0UL, _num_hypernodes), incident_net_prefix_sum);
    });

    ASSERT(pin_prefix_sum.total_sum() == incident_net_prefix_sum.total_sum());
    _num_pins = pin_prefix_sum.total_sum();
    _incident_nets.resize(_num_pins);
    _incidence_array.resize(_num_pins);

    AtomicCounter incident_nets_position(_num_hypernodes,
      parallel::IntegralAtomicWrapper<size_t>(0));
    tbb::parallel_invoke([&] {
      tbb::parallel_for(0UL, _num_hyperedges, [&](const size_t pos) {
        // Setup hyperedges
        Hyperedge& hyperedge = _hyperedges[pos];
        hyperedge.enable();
        hyperedge.setFirstEntry(pin_prefix_sum[pos]);
        hyperedge.setSize(pin_prefix_sum.value(pos));
        if ( hyperedge_weight ) {
          hyperedge.setWeight(hyperedge_weight[pos]);
        }

        const HyperedgeID he = pos;
        size_t incidence_array_pos = hyperedge.firstEntry();
        size_t hash = kEdgeHashSeed;
        for ( const HypernodeID& pin : edge_vector[pos] ) {
          ASSERT(incidence_array_pos < hyperedge.firstInvalidEntry());
          ASSERT(pin < _num_hypernodes);
          // Compute hash of hyperedge
          hash += kahypar::math::hash(pin);
          // Add pin to incidence array
          _incidence_array[incidence_array_pos++] = pin;
          // Add hyperedge he as a incident net to pin
          const size_t incident_nets_pos = incident_net_prefix_sum[pin] +
            incident_nets_position[pin]++;
          ASSERT(incident_nets_pos < incident_net_prefix_sum[pin + 1]);
          _incident_nets[incident_nets_pos] = he;
        }
        hyperedge.hash() = hash;
      });
    }, [&] {
      tbb::parallel_for(0UL, _num_hypernodes, [&](const size_t pos) {
        // Setup hypernodes
        Hypernode& hypernode = _hypernodes[pos];
        hypernode.enable();
        hypernode.setFirstEntry(incident_net_prefix_sum[pos]);
        hypernode.setSize(incident_net_prefix_sum.value(pos));
        if ( hypernode_weight ) {
          hypernode.setWeight(hypernode_weight[pos]);
        }
      });
    });

    // Compute total weight of hypergraph
    _total_weight = tbb::parallel_reduce(
      tbb::blocked_range<HypernodeID>(0UL, _num_hypernodes), 0,
      [&](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
        HypernodeWeight weight = init;
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
          weight += _hypernodes[hn].weight();
        }
        return weight;
      }, std::plus<HypernodeWeight>());
  }

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;

  // ! Hypernodes
  parallel::scalable_vector<Hypernode> _hypernodes;
  // ! Pins of hyperedges
  IncidentNets _incident_nets;
  // ! Hyperedges
  parallel::scalable_vector<Hyperedge> _hyperedges;
  // ! Incident nets of hypernodes
  IncidenceArray _incidence_array;
};

} // namespace ds
} // namespace mt_kahypar
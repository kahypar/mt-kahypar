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


#include <tbb/parallel_for.h>

#include "include/mtkahypartypes.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {
    namespace ds {

// Forward
        class MutableHypergraphFactory;
        template <typename Hypergraph,
                typename ConnectivityInformation>
        class PartitionedHypergraph;

        class MutableHypergraph {

            static constexpr bool enable_heavy_assert = false;

            // During contractions we temporary memcpy all incident nets of a collapsed
            // vertex to consecutive range in a temporary incident nets structure.
            // Afterwards, we sort that range and remove duplicates. However, it turned
            // out that this become a major sequential bottleneck in presence of high
            // degree vertices. Therefore, all vertices with temporary degree greater
            // than this threshold are contracted with a special procedure.
//            static constexpr HyperedgeID HIGH_DEGREE_CONTRACTION_THRESHOLD = ID(500000);

            static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
            static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

//            using AtomicHypernodeID = parallel::IntegralAtomicWrapper<HypernodeID>;
//            using AtomicHypernodeWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;
            using UncontractionFunction = std::function<void (const HypernodeID, const HypernodeID, const HyperedgeID)>;
#define NOOP_BATCH_FUNC [] (const HypernodeID, const HypernodeID, const HyperedgeID) { }

            /**
             * Represents a hypernode of the hypergraph and contains all information
             * associated with a vertex.
             */
            class Hypernode {
            public:
                using IDType = HypernodeID;

                Hypernode() :
                        _incident_nets(),
                        _begin(0),
                        _size(0),
                        _weight(1),
                        _valid(false) { }

                Hypernode(const bool valid) :
                        _incident_nets(),
                        _begin(0),
                        _size(0),
                        _weight(1),
                        _valid(valid) { }

                // Sentinel Constructor
                Hypernode(const size_t begin) :
                        _incident_nets(begin),
                        _begin(begin),
                        _size(0),
                        _weight(1),
                        _valid(false) { }

                Hypernode(const std::vector<HyperedgeID>& incident_nets, const HypernodeWeight weight) :
                        _incident_nets(incident_nets),
                        _begin(0),
                        _size(incident_nets.size()),
                        _weight(weight),
                        _valid(true) { }

                bool isDisabled() const {
                  return !_valid;
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

                // ! returns a pointer to the first element in _incident_nets
                const HyperedgeID* firstEntryPtr() const {
                  return _incident_nets.data() + firstEntry();
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
                  return _weight;
                }

                void setWeight(HyperedgeWeight weight) {
                  ASSERT(!isDisabled());
                  _weight = weight;
                }

                void mark_deleted() {
                  _deleted = true;
                  //TODO maybe also set _valid to false?
                }

                bool is_deleted() const {
                  return _deleted;
                }

                // ! Returns a range to loop over the incident nets of hypernode u.
                IteratorRange<std::vector<HypernodeID>::const_iterator> incidentEdges(const size_t pos = 0) const {
                  return IteratorRange<IncidentNetsIterator>(
                          _incident_nets.cbegin() + firstEntry() + pos,
                          _incident_nets.cbegin() + firstInvalidEntry());
                }

                // ! Adds a net to the hypernode
                void addNet(const HyperedgeID net) {
                  _incident_nets.insert(_incident_nets.begin() + firstInvalidEntry(), net);
                  ++_size;
                }

                // ! Deletes a net from the hypernode
                void deleteNet(const HyperedgeID net) {
                  _incident_nets.erase(
                      std::remove(_incident_nets.begin() + firstEntry(),
                                  _incident_nets.begin() + firstInvalidEntry(), net),
                      _incident_nets.begin() + firstInvalidEntry());
                  --_size;
                }

                // ! activate net
                void activateNet(const HyperedgeID e) {
                  using std::swap;

                  size_t incident_nets_pos = firstInvalidEntry();
                  for ( ; incident_nets_pos > firstEntry(); --incident_nets_pos ) {
                    if ( _incident_nets[incident_nets_pos - 1] == e ) {
                      break;
                    }
                  }
                  ASSERT(incident_nets_pos >= firstEntry());
                  swap(_incident_nets[incident_nets_pos], _incident_nets[firstInvalidEntry() - 1]);
                  _size = _size + 1;
                }

                // ! deactivate net
                void deactivateNet(const HyperedgeID e) {
                  using std::swap;

                  size_t incident_nets_pos = firstEntry();
                  for ( ; incident_nets_pos < firstInvalidEntry(); ++incident_nets_pos ) {
                    if ( _incident_nets[incident_nets_pos] == e ) {
                      break;
                    }
                  }
                  ASSERT(incident_nets_pos < firstInvalidEntry());
                  swap(_incident_nets[incident_nets_pos], _incident_nets[firstInvalidEntry() - 1]);
                  _size = _size - 1;
                }

            private:
                // ! incident_nets
                std::vector<HyperedgeID> _incident_nets;
                // ! Index of the first element in _incident_nets
                size_t _begin;
                // ! Number of incident nets
                size_t _size;
                // ! Hypernode weight
                HypernodeWeight _weight;
                // ! Flag indicating whether or not the element is active.
                bool _valid;
                // ! Flag indicating whether or not the element is deleted.
                 bool _deleted = false;
            };

            /**
             * Represents a hyperedge of the hypergraph and contains all information
             * associated with a net (except connectivity information).
             */
            class Hyperedge {
            public:
                using IDType = HyperedgeID;

                Hyperedge() :
                        _incidence_array(),
                        _begin(0),
                        _size(0),
                        _weight(1),
                        _valid(false) { }

                // Sentinel Constructor
                Hyperedge(const size_t begin) :
                        _incidence_array(begin),
                        _begin(begin),
                        _size(0),
                        _weight(1),
                        _valid(false) { }

                Hyperedge(const std::vector<HypernodeID>& incidence_array, const HyperedgeWeight weight) :
                        _incidence_array(incidence_array),
                        _begin(0),
                        _size(incidence_array.size()),
                        _weight(weight),
                        _valid(true) { }

                Hyperedge(const parallel::scalable_vector<HypernodeID>& incidence_array,
                          const HyperedgeWeight weight) :
                        _incidence_array(incidence_array.begin(), incidence_array.end()),
                        _begin(0),
                        _size(incidence_array.size()),
                        _weight(weight),
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

                // ! returns a pointer to the first element in _incidence_array
                const HypernodeID* firstEntryPtr() const {
                  return _incidence_array.data() + firstEntry();
                }

                // ! Sets the index of the first element in _incidence_array to begin
//                void setFirstEntry(size_t begin) {
//                  ASSERT(!isDisabled());
//                  _begin = begin;
//                }

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

                bool operator== (const Hyperedge& rhs) const {
                  return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
                }

                bool operator!= (const Hyperedge& rhs) const {
                  return _begin != rhs._begin || _size != rhs._size || _weight != rhs._weight;
                }

                void mark_deleted() {
                  _deleted = true;
                }

                bool is_deleted() const {
                  return _deleted;
                }

                // ! Returns a range to loop over the incident nets of the hyperedge
                IteratorRange<std::vector<HypernodeID>::const_iterator> pins() const {
                  return IteratorRange<IncidenceIterator>(
                          _incidence_array.cbegin() + firstEntry(),
                          _incidence_array.cbegin() + firstInvalidEntry());
                }

                //adds a pin to the hyperedge
                void addPin(const HypernodeID pin) {
                  if (firstInvalidEntry() >= _incidence_array.size()) {
                    _incidence_array.push_back(pin);
                  } else {
                    _incidence_array.insert(_incidence_array.begin() + firstInvalidEntry(), pin);
                  }
                  ++_size;
                }

                // ! Deletes a pin from the hyperedge
                void deletePin(const HypernodeID pin) {
                  _incidence_array.erase(
                      std::remove(_incidence_array.begin() + firstEntry(),
                                  _incidence_array.begin() + firstInvalidEntry(), pin),
                      _incidence_array.begin() + firstInvalidEntry());
                  --_size;
                }

                // ! activate pin
                void activatePin(const HypernodeID pin) {
                  using std::swap;

                  size_t incidence_array_pos = firstInvalidEntry();
                  for ( ; incidence_array_pos > firstEntry(); --incidence_array_pos ) {
                    if ( _incidence_array[incidence_array_pos - 1] == pin ) {
                      break;
                    }
                  }
                  ASSERT(incidence_array_pos >= firstEntry());
                  swap(_incidence_array[incidence_array_pos], _incidence_array[firstInvalidEntry() - 1]);
                  _size = _size + 1;
                }

                // ! deactivate pin
                void deactivatePin(const HypernodeID pin) {
                  using std::swap;

                  size_t incidence_array_pos = firstEntry();
                  for ( ; incidence_array_pos < firstInvalidEntry(); ++incidence_array_pos ) {
                    if ( _incidence_array[incidence_array_pos] == pin ) {
                      break;
                    }
                  }
                  ASSERT(incidence_array_pos < firstInvalidEntry());
                  swap(_incidence_array[incidence_array_pos], _incidence_array[firstInvalidEntry() - 1]);
                  _size = _size - 1;
                }

                // ! map pin to a new pin
                void mapPin(const HypernodeID old_pin, const HypernodeID new_pin) {
                  for (size_t pos = firstEntry(); pos < firstInvalidEntry(); ++pos) {
                    if (_incidence_array[pos] == old_pin) {
                      _incidence_array[pos] = new_pin;
                      return;
                    }
                  }
                  ASSERT(false, "Pin " << old_pin << " not found in hyperedge " << _begin);
                }

            private:
                // ! incidence_array
                std::vector<HypernodeID> _incidence_array;
                // ! Index of the first element in _incidence_array
                size_t _begin;
                // ! Number of pins
                size_t _size;
                // ! hyperedge weight
                HyperedgeWeight _weight;
                // ! Flag indicating whether or not the element is active.
                bool _valid;
                // ! Flag indicating whether or not the element is deleted.
                bool _deleted = false;
            };

            /*!
             * Iterator for HypergraphElements (Hypernodes/Hyperedges)
             *
             * The iterator is used in for-each loops over all hypernodes/hyperedges.
             * In order to support iteration over coarsened hypergraphs, this iterator
             * skips over HypergraphElements marked as invalid.
             * Iterating over the set of vertices \f$V\f$ therefore is linear in the
             * size \f$|V|\f$ of the original hypergraph - even if it has been coarsened
             * to much smaller size. The same also holds true for for-each loops over
             * the set of hyperedges.
             *
             * In order to be as generic as possible, the iterator does not expose the
             * internal Hypernode/Hyperedge representations. Instead only handles to
             * the respective elements are returned, i.e. the IDs of the corresponding
             * hypernodes/hyperedges.
             *
             */
            template <typename ElementType>
            class HypergraphElementIterator {
            public:
                using IDType = typename ElementType::IDType;
                using iterator_category = std::forward_iterator_tag;
                using value_type = IDType;
                using reference = IDType&;
                using pointer = const IDType*;
                using difference_type = std::ptrdiff_t;

                /*!
                 * Construct a HypergraphElementIterator
                 * See GenericHypergraph::nodes() or GenericHypergraph::edges() for usage.
                 *
                 * If start_element is invalid, the iterator advances to the first valid
                 * element.
                 *
                 * \param start_element A pointer to the starting position
                 * \param id The index of the element the pointer points to
                 * \param max_id The maximum index allowed
                 */
                HypergraphElementIterator(const ElementType* start_element, IDType id, IDType max_id) :
                        _id(id),
                        _max_id(max_id),
                        _element(start_element) {
                  if (_id != _max_id && _element->isDisabled()) {
                    operator++ ();
                  }
                }

                // ! Returns the id of the element the iterator currently points to.
                IDType operator* () const {
                  return _id;
                }

                // ! Prefix increment. The iterator advances to the next valid element.
                HypergraphElementIterator & operator++ () {
                  ASSERT(_id < _max_id);
                  do {
                    ++_id;
                    ++_element;
                  } while (_id < _max_id && _element->isDisabled());
                  return *this;
                }

                // ! Postfix increment. The iterator advances to the next valid element.
                HypergraphElementIterator operator++ (int) {
                  HypergraphElementIterator copy = *this;
                  operator++ ();
                  return copy;
                }

                bool operator!= (const HypergraphElementIterator& rhs) {
                  return _id != rhs._id;
                }

                bool operator== (const HypergraphElementIterator& rhs) {
                  return _id == rhs._id;
                }

            private:
                // Handle to the HypergraphElement the iterator currently points to
                IDType _id = 0;
                // Maximum allowed index
                IDType _max_id = 0;
                // HypergraphElement the iterator currently points to
                const ElementType* _element = nullptr;
            };


            //TODO did this break something?
//            static_assert(std::is_trivially_copyable<Hypernode>::value, "Hypernode is not trivially copyable");
//            static_assert(std::is_trivially_copyable<Hyperedge>::value, "Hyperedge is not trivially copyable");

            using IncidenceArray = std::vector<HypernodeID>;
            using IncidentNets = std::vector<HyperedgeID>;

            // ! Contains buffers that are needed during multilevel contractions.
            // ! Struct is allocated on top level hypergraph and passed to each contracted
            // ! hypergraph such that memory can be reused in consecutive contractions.
            struct TmpContractionBuffer {
                explicit TmpContractionBuffer(const HypernodeID num_hypernodes,
                                              const HyperedgeID num_hyperedges,
                                              const HyperedgeID num_pins [[maybe_unused]]) {
                  tbb::parallel_invoke([&] {
                      mapping.resize(num_hypernodes);
                  }, [&] {
                      tmp_hypernodes.resize(num_hypernodes);
                  }, [&] {
                      tmp_num_incident_nets.resize(num_hypernodes);
                  }, [&] {
                      hn_weights.resize(num_hypernodes);
                  }, [&] {
                      tmp_hyperedges.resize(num_hyperedges);
                  }, [&] {
                      he_sizes.resize(num_hyperedges);
                  }, [&] {
                      valid_hyperedges.resize(num_hyperedges);
                  });
                }

                std::vector<size_t> mapping;
                std::vector<Hypernode> tmp_hypernodes;
                std::vector<parallel::IntegralAtomicWrapper<size_t>> tmp_num_incident_nets;
                std::vector<parallel::IntegralAtomicWrapper<HypernodeWeight>> hn_weights;
                std::vector<Hyperedge> tmp_hyperedges;
                std::vector<size_t> he_sizes;
                std::vector<size_t> valid_hyperedges;
            };

        public:
            static constexpr bool is_graph = false;
            static constexpr bool is_static_hypergraph = true;
            static constexpr bool is_partitioned = false;
            static constexpr size_t SIZE_OF_HYPERNODE = sizeof(Hypernode);
            static constexpr size_t SIZE_OF_HYPEREDGE = sizeof(Hyperedge);
            static constexpr mt_kahypar_hypergraph_type_t TYPE = MUTABLE_HYPERGRAPH;

            // ! Factory
            using Factory = MutableHypergraphFactory;
            // ! Iterator to iterate over the hypernodes
            using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
            // ! Iterator to iterate over the hyperedges
            using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
            // ! Iterator to iterate over the pins of a hyperedge
            using IncidenceIterator = typename IncidenceArray::const_iterator;
            // ! Iterator to iterate over the incident nets of a hypernode
            using IncidentNetsIterator = typename IncidentNets::const_iterator;

            struct ParallelHyperedge {
                HyperedgeID removed_hyperedge;
                HyperedgeID representative;
            };

            explicit MutableHypergraph() :
                    _num_hypernodes(0),
                    _num_removed_hypernodes(0),
                    _removed_degree_zero_hn_weight(0),
                    _num_hyperedges(0),
                    _num_removed_hyperedges(0),
                    _max_edge_size(0),
                    _num_pins(0),
                    _total_degree(0),
                    _total_weight(0),
                    _hypernodes(),
                    _hyperedges(),
                    _community_ids(0),
                    _fixed_vertices(),
                    _tmp_contraction_buffer(nullptr) { }

            MutableHypergraph(const MutableHypergraph&) = delete;
            MutableHypergraph & operator= (const MutableHypergraph &) = delete;

            MutableHypergraph(MutableHypergraph&& other) :
                    _num_hypernodes(other._num_hypernodes),
                    _num_removed_hypernodes(other._num_removed_hypernodes),
                    _removed_degree_zero_hn_weight(other._removed_degree_zero_hn_weight),
                    _num_hyperedges(other._num_hyperedges),
                    _num_removed_hyperedges(other._num_removed_hyperedges),
                    _max_edge_size(other._max_edge_size),
                    _num_pins(other._num_pins),
                    _total_degree(other._total_degree),
                    _total_weight(other._total_weight),
                    _hypernodes(std::move(other._hypernodes)),
                    _hyperedges(std::move(other._hyperedges)),
                    _community_ids(std::move(other._community_ids)),
                    _fixed_vertices(std::move(other._fixed_vertices)),
                    _tmp_contraction_buffer(std::move(other._tmp_contraction_buffer)) {
              _fixed_vertices.setHypergraph(this);
              other._tmp_contraction_buffer = nullptr;
            }

            MutableHypergraph & operator= (MutableHypergraph&& other) {
              _num_hypernodes = other._num_hypernodes;
              _num_removed_hypernodes = other._num_removed_hypernodes;
              _removed_degree_zero_hn_weight = other._removed_degree_zero_hn_weight;
              _num_hyperedges = other._num_hyperedges;
              _num_removed_hyperedges = other._num_removed_hyperedges;
              _max_edge_size = other._max_edge_size;
              _num_pins = other._num_pins;
              _total_degree = other._total_degree;
              _total_weight = other._total_weight;
              _hypernodes = std::move(other._hypernodes);
              _hyperedges = std::move(other._hyperedges);
              _community_ids = std::move(other._community_ids);
              _fixed_vertices = std::move(other._fixed_vertices);
              _fixed_vertices.setHypergraph(this);
              _tmp_contraction_buffer = std::move(other._tmp_contraction_buffer);
              other._tmp_contraction_buffer = nullptr;
              return *this;
            }

            ~MutableHypergraph() {
              if ( _tmp_contraction_buffer ) {
                delete(_tmp_contraction_buffer);
                _tmp_contraction_buffer = nullptr;
              }
              freeInternalData();
            }

            // ####################### General Hypergraph Stats #######################

            // ! Initial number of hypernodes
            HypernodeID initialNumNodes() const {
              return _num_hypernodes;
            }

            // ! Number of removed hypernodes
            HypernodeID numRemovedHypernodes() const {
              return _num_removed_hypernodes;
            }

            // ! Weight of removed degree zero vertics
            HypernodeWeight weightOfRemovedDegreeZeroVertices() const {
              return _removed_degree_zero_hn_weight;
            }

            // ! Initial number of hyperedges
            HyperedgeID initialNumEdges() const {
              return _num_hyperedges;
            }

            // ! Number of removed hyperedges
            HyperedgeID numRemovedHyperedges() const {
              return _num_removed_hyperedges;
            }

            // ! Set the number of removed hyperedges
            void setNumRemovedHyperedges(const HyperedgeID num_removed_hyperedges) {
              _num_removed_hyperedges = num_removed_hyperedges;
            }

            // ! Initial number of pins
            HypernodeID initialNumPins() const {
              return _num_pins;
            }

            // ! Initial sum of the degree of all vertices
            HypernodeID initialTotalVertexDegree() const {
              return _total_degree;
            }

            // ! Total weight of hypergraph
            HypernodeWeight totalWeight() const {
              return _total_weight;
            }

            // ! Computes the total node weight of the hypergraph
            void computeAndSetTotalNodeWeight(parallel_tag_t);

            // ####################### Iterators #######################

            // ! Iterates in parallel over all active nodes and calls function f
            // ! for each vertex
            template<typename F>
            void doParallelForAllNodes(const F& f) const {
              tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
                  if ( nodeIsEnabled(hn) ) {
                    f(hn);
                  }
              });
            }

            // ! Iterates in parallel over all active edges and calls function f
            // ! for each net
            template<typename F>
            void doParallelForAllEdges(const F& f) const {
              tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
                  if ( edgeIsEnabled(he) ) {
                    f(he);
                  }
              });
            }

            // ! Returns a range of the active nodes of the hypergraph
            IteratorRange<HypernodeIterator> nodes() const {
              return IteratorRange<HypernodeIterator>(
                      HypernodeIterator(_hypernodes.data(), ID(0), _num_hypernodes),
                      HypernodeIterator(_hypernodes.data() + _num_hypernodes, _num_hypernodes, _num_hypernodes));
            }

            // ! Returns a range of the active edges of the hypergraph
            IteratorRange<HyperedgeIterator> edges() const {
              return IteratorRange<HyperedgeIterator>(
                      HyperedgeIterator(_hyperedges.data(), ID(0), _num_hyperedges),
                      HyperedgeIterator(_hyperedges.data() + _num_hyperedges, _num_hyperedges, _num_hyperedges));
            }

            // ! Returns a range to loop over the incident nets of hypernode u.
            IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
              ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
              const Hypernode& hn = hypernode(u);
              return hn.incidentEdges();
            }

            // ! Returns a range to loop over the pins of hyperedge e.
            IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
              ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
              const Hyperedge& he = hyperedge(e);
              return he.pins();
            }

            // ####################### Hypernode Information #######################

            // ! Weight of a vertex
            HypernodeWeight nodeWeight(const HypernodeID u) const {
              return hypernode(u).weight();
            }

            // ! Sets the weight of a vertex
            void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
              ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
              return hypernode(u).setWeight(weight);
            }

            // ! Degree of a hypernode
            HyperedgeID nodeDegree(const HypernodeID u) const {
              ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
              return hypernode(u).size();
            }

            // ! Returns, whether a hypernode is enabled or not
            bool nodeIsEnabled(const HypernodeID u) const {
              return !hypernode(u).isDisabled();
            }

            // ! Enables a hypernode (must be disabled before)
            void enableHypernode(const HypernodeID u) {
              hypernode(u).enable();
            }

            // ! Disables a hypernode (must be enabled before)
            void disableHypernode(const HypernodeID u) {
              hypernode(u).disable();
            }

            // ! Removes a hypernode (must be enabled before)
            void removeHypernode(const HypernodeID u) {
              hypernode(u).disable();
              ++_num_removed_hypernodes;
            }

            // ! Removes a degree zero hypernode
            void removeDegreeZeroHypernode(const HypernodeID u) {
              ASSERT(nodeDegree(u) == 0);
              _removed_degree_zero_hn_weight += nodeWeight(u);
              removeHypernode(u);
            }

            // ! Restores a degree zero hypernode
            void restoreDegreeZeroHypernode(const HypernodeID u) {
              hypernode(u).enable();
              ASSERT(nodeDegree(u) == 0);
              _removed_degree_zero_hn_weight -= nodeWeight(u);
            }

            // ####################### Hyperedge Information #######################

            // ! Weight of a hyperedge
            HypernodeWeight edgeWeight(const HyperedgeID e) const {
              ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
              return hyperedge(e).weight();
            }

            // ! Sets the weight of a hyperedge
            void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
              ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
              return hyperedge(e).setWeight(weight);
            }

            // ! Number of pins of a hyperedge
            HypernodeID edgeSize(const HyperedgeID e) const {
              ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
              return hyperedge(e).size();
            }

            // ! Maximum size of a hyperedge
            HypernodeID maxEdgeSize() const {
              return _max_edge_size;
            }

            // ! Returns, whether a hyperedge is enabled or not
            bool edgeIsEnabled(const HyperedgeID e) const {
              return !hyperedge(e).isDisabled();
            }

            // ! Enables a hyperedge (must be disabled before)
            void enableHyperedge(const HyperedgeID e) {
              hyperedge(e).enable();
            }

            // ! Disabled a hyperedge (must be enabled before)
            void disableHyperedge(const HyperedgeID e) {
              hyperedge(e).disable();
            }

            // ! Community id which hypernode u is assigned to
            PartitionID communityID(const HypernodeID u) const {
              return _community_ids[u];
            }

            // ! Assign a community to a hypernode
            void setCommunityID(const HypernodeID u, const PartitionID community_id) {
              _community_ids[u] = community_id;
            }

            // ####################### Fixed Vertex Support #######################

            void addFixedVertexSupport(FixedVertexSupport<MutableHypergraph>&& fixed_vertices) {
              _fixed_vertices = std::move(fixed_vertices);
              _fixed_vertices.setHypergraph(this);
            }

            bool hasFixedVertices() const {
              return _fixed_vertices.hasFixedVertices();
            }

            HypernodeWeight totalFixedVertexWeight() const {
              return _fixed_vertices.totalFixedVertexWeight();
            }

            HypernodeWeight fixedVertexBlockWeight(const PartitionID block) const {
              return _fixed_vertices.fixedVertexBlockWeight(block);
            }

            bool isFixed(const HypernodeID hn) const {
              return _fixed_vertices.isFixed(hn);
            }

            PartitionID fixedVertexBlock(const HypernodeID hn) const {
              return _fixed_vertices.fixedVertexBlock(hn);
            }

            void setMaxFixedVertexBlockWeight(const std::vector<HypernodeWeight> max_block_weights) {
              _fixed_vertices.setMaxBlockWeight(max_block_weights);
            }

            const FixedVertexSupport<MutableHypergraph>& fixedVertexSupport() const {
              return _fixed_vertices;
            }

            FixedVertexSupport<MutableHypergraph> copyOfFixedVertexSupport() const {
              return _fixed_vertices.copy();
            }

            // ####################### Contract / Uncontract #######################

            /*!
             * Contracts a given community structure. All vertices with the same label
             * are collapsed into the same vertex. The resulting single-pin and parallel
             * hyperedges are removed from the contracted graph. The function returns
             * the contracted hypergraph and a mapping which specifies a mapping from
             * community label (given in 'communities') to a vertex in the coarse hypergraph.
             *
             * \param communities Community structure that should be contracted
             */
            MutableHypergraph contract(parallel::scalable_vector<HypernodeID>& communities, bool deterministic = false);

            bool registerContraction(const HypernodeID, const HypernodeID) {
              throw UnsupportedOperationException(
                      "registerContraction(u, v) is not supported in static hypergraph");
              return false;
            }

            size_t contract(const HypernodeID,
                            const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max()) {
              unused(max_node_weight);
              throw UnsupportedOperationException(
                      "contract(v, max_node_weight) is not supported in static hypergraph");
              return 0;
            }

            void uncontract(const Batch&,
                            const UncontractionFunction& case_one_func = NOOP_BATCH_FUNC,
                            const UncontractionFunction& case_two_func = NOOP_BATCH_FUNC) {
              unused(case_one_func);
              unused(case_two_func);
              throw UnsupportedOperationException(
                      "uncontract(batch) is not supported in static hypergraph");
            }

            VersionedBatchVector createBatchUncontractionHierarchy(const size_t) {
              throw UnsupportedOperationException(
                      "createBatchUncontractionHierarchy(batch_size) is not supported in static hypergraph");
              return { };
            }

            // ####################### Remove / Restore Hyperedges #######################

            /*!
            * Removes a hyperedge from the hypergraph. This includes the removal of he from all
            * of its pins and to disable the hyperedge.
            *
            * NOTE, this function is not thread-safe and should only be called in a single-threaded
            * setting.
            */
            void removeEdge(const HyperedgeID he) {
              ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
              for ( const HypernodeID& pin : pins(he) ) {
                removeIncidentEdgeFromHypernode(he, pin);
              }
              ++_num_removed_hyperedges;
              disableHyperedge(he);
            }

            /*!
            * Removes a hyperedge from the hypergraph. This includes the removal of he from all
            * of its pins and to disable the hyperedge. Noze, in contrast to removeEdge, this function
            * removes hyperedge from all its pins in parallel.
            *
            * NOTE, this function is not thread-safe and should only be called in a single-threaded
            * setting.
            */
            void removeLargeEdge(const HyperedgeID he) {
              ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
//              const size_t incidence_array_start = hyperedge(he).firstEntry();
//              const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
//              tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
//                  const HypernodeID pin = _incidence_array[pos];
//                  removeIncidentEdgeFromHypernode(he, pin);
//              });
              disableHyperedge(he);
            }

            /*!
             * Restores a large hyperedge previously removed from the hypergraph.
             */
            void restoreLargeEdge(const HyperedgeID& he) {
              ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "is enabled");
              enableHyperedge(he);
//              const size_t incidence_array_start = hyperedge(he).firstEntry();
//              const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
//              tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
//                  const HypernodeID pin = _incidence_array[pos];
//                  insertIncidentEdgeToHypernode(he, pin);
//              });
            }

            parallel::scalable_vector<ParallelHyperedge> removeSinglePinAndParallelHyperedges() {
              throw UnsupportedOperationException(
                      "removeSinglePinAndParallelHyperedges() is not supported in static hypergraph");
              return { };
            }

            void restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>&) {
              throw UnsupportedOperationException(
                      "restoreSinglePinAndParallelNets(hes_to_restore) is not supported in static hypergraph");
            }

            // ####################### Initialization / Reset Functions #######################

            // ! Reset internal community information
            void copyCommunityIDs(const parallel::scalable_vector<PartitionID>& community_ids) {
              ASSERT(community_ids.size() == UI64(_num_hypernodes));
              doParallelForAllNodes([&](const HypernodeID& hn) {
                  _community_ids[hn] = community_ids[hn];
              });
            }

            void setCommunityIDs(ds::Clustering&& communities) {
              ASSERT(communities.size() == initialNumNodes());
              _community_ids = std::move(communities);
            }

            // ! Copy static hypergraph in parallel
            MutableHypergraph copy(parallel_tag_t) const;

            // ! Copy static hypergraph sequential
            MutableHypergraph copy() const;

            // ! Convert to static hypergraph
            StaticHypergraph toStaticHypergraph(std::vector<HypernodeID>* old_to_new_hn,
                                                std::vector<HyperedgeID>* old_to_new_he, std::vector<HypernodeID>* deleted_hn, std::vector<HyperedgeID>* deleted_he);

            // ! Convert StaticHypergraph to MutableHypergraph
            MutableHypergraph fromStaticHypergraph(const StaticHypergraph& static_hg, std::vector<HypernodeID>* old_to_new_hn,
                                                   std::vector<HyperedgeID>* old_to_new_he, std::vector<HypernodeID>* deleted_hn, std::vector<HyperedgeID>* deleted_he);

            // ! Reset internal data structure
            void reset() { }

            // ! Free internal data in parallel
            void freeInternalData() {
              if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
                freeTmpContractionBuffer();
              }
              _num_hypernodes = 0;
              _num_hyperedges = 0;
            }

            void freeTmpContractionBuffer() {
              if ( _tmp_contraction_buffer ) {
                delete(_tmp_contraction_buffer);
                _tmp_contraction_buffer = nullptr;
              }
            }

            void memoryConsumption(utils::MemoryTreeNode* parent) const;

            // ! Only for testing
            bool verifyIncidenceArrayAndIncidentNets() {
              throw UnsupportedOperationException(
                      "verifyIncidenceArrayAndIncidentNets() not supported in static hypergraph");
              return false;
            }

            // ######################## Mutations ########################

            /*!
             * Adds a new hypernode to the hypergraph. The hypernode is initialized with
             * the given weight and is assigned the next available hypernode ID.
             *
             * \param hypernode The hypernode to be added
             */
            HypernodeID addHypernode(const std::vector<HyperedgeID>& incident_nets, const HypernodeWeight weight = 1) {
              Hypernode hypernode({}, weight);
              //insert new hypernode in front of sentinel
              _hypernodes.insert(_hypernodes.end() - 1, hypernode);
              _community_ids.push_back(0);
              ++_num_hypernodes;
              _total_weight += hypernode.weight();
              for (const HyperedgeID& e : incident_nets) {
                addPin(e, _num_hypernodes - 1);
              }
              return _num_hypernodes - 1;
            }

            /*!
             * Deletes a hypernode from the hypergraph. This includes the removal of all its adjacancies
              * from the hyperedges and the incident nets.
              * \param u The hypernode to be deleted
              */
            void deleteHypernode(const HypernodeID u) {
              ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
              ASSERT(!hypernode(u).is_deleted(), "Hypernode" << u << "is deleted");
              _num_hypernodes--;
              _total_weight -= hypernode(u).weight();
              for (const HyperedgeID& e : incident_nets_of(u)) {
                deletePin(e, u);
              }
              hypernode(u).mark_deleted();
            }

            /*!
             * Adds a new hyperedge to the hypergraph. The hyperedge is initialized with
             * the given weight and is assigned the next available hyperedge ID.
             *
             * \param hyperedge The hyperedge to be added
             */
            HyperedgeID addHyperedge(const std::vector<HypernodeID>& pins,
                                     const HyperedgeWeight weight = 1) {
              Hyperedge hyperedge(pins, weight);
              //insert new hyperedge in front of sentinel
              _hyperedges.insert(_hyperedges.end() - 1, hyperedge);
              ++_num_hyperedges;
              _max_edge_size = std::max(static_cast<size_t>(_max_edge_size), hyperedge.size());
              for (const HypernodeID& pin : pins) {
                addPin(_num_hyperedges - 1, pin);
              }
              return _num_hyperedges - 1;
            }

            /*!
             * Removes a hyperedge from the hypergraph
             * \param he The hyperedge to be removed
             */
            void deleteHyperedge(const HyperedgeID he) {
              ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
              ASSERT(!hyperedge(he).is_deleted(), "Hyperedge" << he << "is deleted");
              _num_hyperedges--;
              for (const HypernodeID& pin : hyperedge(he).pins()) {
                deletePin(he, pin);
              }
              hyperedge(he).mark_deleted();
            }

            /*!
             * Connects a hyperedge to a hypernode by adding the pin to the incidence array
             * and the hyperedge to the incident nets of the hypernode.
             * This function assumes that the hyperedge  and the hypernode are enabled.
             * \param he The hyperedge to which the pin should be added
             * \param pin The pin to be added
             */
            void addPin(const HyperedgeID he, const HypernodeID pin) {
              ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
              ASSERT(!hypernode(pin).isDisabled(), "Hypernode" << pin << "is disabled");
              ASSERT(!hyperedge(he).is_deleted(), "Hyperedge" << he << "is deleted");
              ASSERT(!hypernode(pin).is_deleted(), "Hypernode" << pin << "is deleted");
              // Add pin to the hyperedge
              hyperedge(he).addPin(pin);
              ++_num_pins;
              _max_edge_size = std::max(static_cast<size_t>(_max_edge_size), hyperedge(he).size());
              _total_degree += 1;
              hypernode(pin).addNet(he);
            }
            /*!
             * Deletes a pin from a hyperedge and removes the hyperedge from the incident nets of the hypernode.
             * This function assumes that the hyperedge and the hypernode are enabled.
             * \param he The hyperedge from which the pin should be deleted
             * \param pin The pin to be deleted
             */
            void deletePin(const HyperedgeID he, const HypernodeID pin) {
              ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
              ASSERT(!hypernode(pin).isDisabled(), "Hypernode" << pin << "is disabled");
              ASSERT(!hyperedge(he).is_deleted(), "Hyperedge" << he << "is deleted");
              ASSERT(!hypernode(pin).is_deleted(), "Hypernode" << pin << "is deleted");
              --_num_pins;
              _total_degree -= 1;
              hypernode(pin).deleteNet(he);
              hyperedge(he).deletePin(pin);
            }

        private:
            friend class MutableHypergraphFactory;
            template<typename Hypergraph>
            friend class CommunitySupport;
            template <typename Hypergraph,
                    typename ConnectivityInformation>
            friend class PartitionedHypergraph;

            // ####################### Hypernode Information #######################

            // ! Accessor for hypernode-related information
            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
              ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
              return _hypernodes[u];
            }

            // ! To avoid code duplication we implement non-const version in terms of const version
            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
              return const_cast<Hypernode&>(static_cast<const MutableHypergraph&>(*this).hypernode(u));
            }

            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE IteratorRange<IncidentNetsIterator> incident_nets_of(const HypernodeID u,
                                                                                                    const size_t pos = 0) const {
              ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
              const Hypernode& hn = hypernode(u);
              return hn.incidentEdges(pos);
            }

            // ####################### Hyperedge Information #######################

            // ! Accessor for hyperedge-related information
            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
              ASSERT(e <= _num_hyperedges, "Hyperedge" << e << "does not exist");
              return _hyperedges[e];
            }

            // ! To avoid code duplication we implement non-const version in terms of const version
            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
              return const_cast<Hyperedge&>(static_cast<const MutableHypergraph&>(*this).hyperedge(e));
            }

            // ####################### Remove / Restore Hyperedges #######################

            // ! Removes hyperedge e from the incident nets of vertex hn
            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID e,
                                                                                    const HypernodeID u) {
              hypernode(u).deactivateNet(e);
            }

            // ! Inserts hyperedge he to incident nets array of vertex hn
            MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernode(const HyperedgeID e,
                                                                                  const HypernodeID u) {
              hypernode(u).activateNet(e);
            }

            // ! Allocate the temporary contraction buffer
            void allocateTmpContractionBuffer() {
              if ( !_tmp_contraction_buffer ) {
                _tmp_contraction_buffer = new TmpContractionBuffer(
                        _num_hypernodes, _num_hyperedges, _num_pins);
              }
            }

            // ! Number of hypernodes
            HypernodeID _num_hypernodes;
            // ! Number of removed hypernodes
            HypernodeID _num_removed_hypernodes;
            // ! Number of removed degree zero hypernodes
            HypernodeWeight _removed_degree_zero_hn_weight;
            // ! Number of hyperedges
            HyperedgeID _num_hyperedges;
            // ! Number of removed hyperedges
            HyperedgeID _num_removed_hyperedges;
            // ! Maximum size of a hyperedge
            HypernodeID _max_edge_size;
            // ! Number of pins
            HypernodeID _num_pins;
            // ! Total degree of all vertices
            HypernodeID _total_degree;
            // ! Total weight of hypergraph
            HypernodeWeight _total_weight;

            // ! Hypernodes
            std::vector<Hypernode> _hypernodes;
            // ! Hyperedges
            std::vector<Hyperedge> _hyperedges;

            // ! Communities
            ds::Clustering _community_ids;

            // ! Fixed Vertex Support
            FixedVertexSupport<MutableHypergraph> _fixed_vertices;

            // ! Data that is reused throughout the multilevel hierarchy
            // ! to contract the hypergraph and to prevent expensive allocations
            TmpContractionBuffer* _tmp_contraction_buffer;
        };

    } // namespace ds
} // namespace mt_kahypar

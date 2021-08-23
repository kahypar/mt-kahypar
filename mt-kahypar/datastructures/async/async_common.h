//
// Created by mlaupichler on 28.04.21.
//

#ifndef KAHYPAR_ASYNC_COMMON_H
#define KAHYPAR_ASYNC_COMMON_H

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "tbb/blocked_range.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/datastructures/array.h"

#include "boost/iterator/iterator_adaptor.hpp"

namespace mt_kahypar::ds {

    using ContractionGroupID = HypernodeID;
    using ContractionGroupIDIterator = parallel::scalable_vector<ContractionGroupID>::const_iterator;
    using ContractionGroupIDIteratorRange = IteratorRange<ContractionGroupIDIterator>;
    using BlockedGroupIDIterator = tbb::blocked_range<ContractionGroupID>;
    using Contraction = Memento;
    using ContractionIterator = std::vector<Contraction>::const_iterator;

    static constexpr ContractionGroupID invalidGroupID = std::numeric_limits<ContractionGroupID>::max();
    static constexpr HypernodeID invalidDepth = std::numeric_limits<HypernodeID>::max();

    /// Represents a group of (un-)contractions that have been contracted simultaneously with the same representative.
    /// This means the contractions in a group have to be uncontracted simultaneously as well in the sense that no local
    /// refinement is performed in between uncontractions of the same group as those states may be inconsistent.
    /// Groups are expected to be small (constant against the number of uncontractions in a graph version).
    class ContractionGroup {

    private:

        std::vector<Contraction> _contractions;
        const HypernodeID _representative;

        /// Check that group is not empty and that all contractions have the same representative
        bool sanityCheck() const {
            if (empty()) return false;
            auto repr = _contractions.front().u;
            return std::all_of(_contractions.begin(),_contractions.end(),[repr](Contraction c){return c.u == repr;});
        };

//        static std::vector<Contraction> sortInitializerList(std::initializer_list<Contraction>& init) {
//            std::vector<Contraction> init_vec(init);
//            std::sort(init_vec.begin(),init_vec.end());
//            return init_vec;
//        }
//
//        static std::vector<Contraction>& sortInitializerVector(std::vector<Contraction>& init) {
//            std::sort(init.begin(), init.end());
//            return init;
//        }

        void sortContractions() {
            std::sort(_contractions.begin(), _contractions.end());
        }

        HypernodeID extractRepresentative() {
            ASSERT(sanityCheck());
            return _contractions.front().u;
        }

    public:
        ContractionIterator begin() {return _contractions.begin();};
        ContractionIterator end() {return _contractions.end();};
        auto cbegin() const {return _contractions.begin();};
        auto cend() const {return _contractions.end();};
        auto begin() const {return _contractions.begin();};
        auto end() const {return _contractions.end();};

        bool empty() const {return _contractions.empty();};
        HypernodeID size() const {return _contractions.size();};

        Contraction at(const int i) const {return _contractions.at(i);};

        ContractionGroup(std::initializer_list<Contraction> init) : _contractions(init), _representative(extractRepresentative()) {
            sortContractions();
        };

        explicit ContractionGroup(std::vector<Contraction>&& init) : _contractions(std::move(init)), _representative(extractRepresentative()) {
            sortContractions();
        };

        ContractionGroup(ContractionGroup&& other)  noexcept :
            _contractions(std::move(other._contractions)),
            _representative(other._representative) {}

        ContractionGroup(const ContractionGroup& other) = delete;
        ContractionGroup& operator=(const ContractionGroup& other) = delete;
        ContractionGroup& operator=(ContractionGroup&& other) = delete;

        HypernodeID getRepresentative() const {
            return _representative;
        }

        /// This function is linear in the number of contractions.
        bool contains(const Contraction& contraction) const {
          return std::find(_contractions.begin(), _contractions.end(), contraction) != _contractions.end();
        }

        bool operator==(const ContractionGroup &rhs) const {
          return _contractions == rhs._contractions;
        }

        bool operator!=(const ContractionGroup &rhs) const {
            return !(rhs == *this);
        }

        void debugPrint() const {
            std::cout << "\tGroup has " << size() << " contractions: \n";
            for (auto &c : *this) {
                std::cout << "\t\t(u: " << c.u << ", v: " << c.v << ")\n";
            }
        }

        size_t size_in_bytes() const {
          return _contractions.size() * sizeof(Contraction) + sizeof(HypernodeID);
        }
    };

    /// Adapts an ContractionIterator to an iterator over the contracted nodes in the contractions.
    /// I.e. this iterates over the 'v' value in the contractions that the given iterator defines.
    class ContractionToNodeIDIteratorAdaptor
        : public boost::iterator_adaptor<
                ContractionToNodeIDIteratorAdaptor,
                ContractionIterator,
                const HypernodeID>
        {

        using iterator = ContractionToNodeIDIteratorAdaptor;

        public:
            ContractionToNodeIDIteratorAdaptor() : ContractionToNodeIDIteratorAdaptor::iterator_adaptor_ () {}
            explicit ContractionToNodeIDIteratorAdaptor(const iterator::iterator_adaptor_::base_type& t) : ContractionToNodeIDIteratorAdaptor::iterator_adaptor_(t) {}

        private:
            friend class boost::iterator_core_access;

            // ! Adapting happens here: Dereference by redirecting to the contracted vertex
            const HypernodeID& dereference() const {return this->base()->v;}
    };

    /// An iterator through all HypernodeIDs in a ContractionGroup, i.e. the representative and the contracted vertices.
    class GroupNodeIDIterator : public std::iterator<
                                        std::forward_iterator_tag, // category
                                        HypernodeID,    // value type
                                        std::ptrdiff_t, // difference type
                                        const HypernodeID*, //pointer type
                                        HypernodeID // reference type
                                        > {
    public:

        // Factory method that gets a GroupNodeIDIterator at the beginning of a given group.
        static GroupNodeIDIterator getAtBegin(const ContractionGroup& group) {
            return GroupNodeIDIterator(&group, false);
        }

        // Factory method that gets a GroupNodeIDIterator at the end of a given group.
        static GroupNodeIDIterator getAtEnd(const ContractionGroup& group) {
            return GroupNodeIDIterator(&group, true);
        }

        // Prefix increment
        GroupNodeIDIterator& operator++() {
            if (_atRepresentative) {
                _atRepresentative = false;
            } else {
                ++_contracted_it;
            }
            return *this;
        }

        // Postfix increment
        GroupNodeIDIterator operator++(int) {
            GroupNodeIDIterator retval = *this;
            ++(*this);
            return retval;
        }

        // Equality operator. Iterators are expected to reference the same group in memory.
        bool operator==(const GroupNodeIDIterator& other) const {
            if (_atRepresentative != other._atRepresentative) return false;
            if (_atRepresentative && other._atRepresentative) {
                if (_group != other._group) return false;
                return true;
            }
            return _contracted_it == other._contracted_it;
        }

        // Inequality operator
        bool operator!=(const GroupNodeIDIterator& other) const {
            return !(*this == other);
        }

        // Dereference operator
        HypernodeID operator*() const {
            if (_atRepresentative) return _group->getRepresentative();
            return *_contracted_it;
        }

    private:

        GroupNodeIDIterator(const ContractionGroup* group, bool atEnd) : _group(group), _atRepresentative(!atEnd) {
            if (atEnd) {
                _contracted_it = ContractionToNodeIDIteratorAdaptor(_group->end());
            } else {
                _contracted_it = ContractionToNodeIDIteratorAdaptor(_group->begin());
            }
        }

        // ! Underlying contraction group
        const ContractionGroup* _group;

        bool _atRepresentative;
        ContractionToNodeIDIteratorAdaptor _contracted_it;
    };

    /// An iterator through a pin snapshot taken for asynchronous uncoarsening. Stitches together a iterator through
    /// stable pins given by a range in the incidence array and an iterator through a vector that contains a copy of
    /// volatile pins.
    class PinSnapshotIterator : public std::iterator<
            std::forward_iterator_tag, // category
            HypernodeID,    // value type
            std::ptrdiff_t, // difference type
            const HypernodeID*, //pointer type
            HypernodeID // reference type
    > {
    private:
        using IncidenceIterator = typename Array<HypernodeID>::iterator;
        using ConstIncidenceIterator = typename Array<HypernodeID>::const_iterator;

    public:

        // Factory method that gets a PinSnapshotIterator range from given stable and volatile ranges
        static IteratorRange<PinSnapshotIterator> stitchPinIterators(IteratorRange<ConstIncidenceIterator> stable_pins_range,
                                                      IteratorRange<HypernodeID*> volatile_pins_range) {
          return IteratorRange<PinSnapshotIterator>(
                  PinSnapshotIterator(stable_pins_range.begin(),stable_pins_range.end(),volatile_pins_range.begin()),
                  PinSnapshotIterator(volatile_pins_range.end(), stable_pins_range.end(), volatile_pins_range.begin()));
        }

        // Prefix increment
        PinSnapshotIterator& operator++() {
          if (_it_in_stable_range) {
            ASSERT(_stable_it != _end_of_stable);
            ++_stable_it;
            if (_stable_it == _end_of_stable) {
              _it_in_stable_range = false;
              ASSERT(_volatile_it == _begin_of_volatile);
            }
          } else {
            ++_volatile_it;
          }
          return *this;
        }

        // Postfix increment
        PinSnapshotIterator operator++(int) {
          PinSnapshotIterator retval = *this;
          ++(*this);
          return retval;
        }

        // Equality operator.
        bool operator==(const PinSnapshotIterator& other) const {
          if (_end_of_stable != other._end_of_stable) return false;
          if (_begin_of_volatile != other._begin_of_volatile) return false;
          return (_stable_it == other._stable_it && _volatile_it == other._volatile_it);
        }

        // Inequality operator
        bool operator!=(const PinSnapshotIterator& other) const {
          return !(*this == other);
        }

        // Dereference operator
        HypernodeID operator*() const {
          if (_it_in_stable_range) {
            ASSERT(_stable_it != _end_of_stable);
            return *_stable_it;
          } else {
            return *_volatile_it;
          }
        }

    private:

        // Constructs a iterator at begin_of_stable
        PinSnapshotIterator(const ConstIncidenceIterator& begin_of_stable,
                            const ConstIncidenceIterator& end_of_stable,
                            HypernodeID* const begin_of_volatile) :
                            _it_in_stable_range(begin_of_stable != end_of_stable),
                            _stable_it(begin_of_stable),
                            _volatile_it(begin_of_volatile),
                            _end_of_stable(end_of_stable),
                            _begin_of_volatile(begin_of_volatile) {}

        // Constructs a iterator at end_of_volatile
        PinSnapshotIterator(HypernodeID* const end_of_volatile,
                            const ConstIncidenceIterator& end_of_stable,
                            HypernodeID* const begin_of_volatile) :
                _it_in_stable_range(false),
                _stable_it(end_of_stable),
                _volatile_it(end_of_volatile),
                _end_of_stable(end_of_stable),
                _begin_of_volatile(begin_of_volatile) {}

        bool _it_in_stable_range;
        IncidenceIterator _stable_it;
        HypernodeID* _volatile_it;

        const IncidenceIterator _end_of_stable;
        HypernodeID* const _begin_of_volatile;
    };

}

#endif //KAHYPAR_ASYNC_COMMON_H

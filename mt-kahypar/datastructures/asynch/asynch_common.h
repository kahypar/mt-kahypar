//
// Created by mlaupichler on 28.04.21.
//

#ifndef KAHYPAR_ASYNCH_COMMON_H
#define KAHYPAR_ASYNCH_COMMON_H

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "tbb/blocked_range.h"
#include "mt-kahypar/utils/range.h"

#include "boost/iterator/iterator_adaptor.hpp"

namespace mt_kahypar::ds {

    using ContractionGroupID = HypernodeID;
    using ContractionGroupIDIteratorRange = IteratorRange<parallel::scalable_vector<ContractionGroupID>::const_iterator>;
    using BlockedGroupIDIterator = tbb::blocked_range<ContractionGroupID>;
    using Contraction = Memento;
    using ContractionIterator = std::vector<Contraction>::const_iterator;


    using IGroupLockManager = ILockManager<HypernodeID,ContractionGroupID>;
    using MockGroupLockManager = MockLockManager<HypernodeID, ContractionGroupID>;

    static constexpr ContractionGroupID invalidGroupID = std::numeric_limits<ContractionGroupID>::max();


    /// Represents a group of (un-)contractions that have been contracted simultaneously with the same representative.
    /// This means the contractions in a group have to be uncontracted simultaneously as well in the sense that no local
    /// refinement is performed in between uncontractions of the same group as those states may be inconsistent.
    /// Groups are expected to be small (constant against the number of uncontractions in a graph version).
    class ContractionGroup {

    private:

        std::vector<Contraction> _contractions;
        HypernodeID _representative;

        /// Check that group is not empty and that all contractions have the same representative
        bool sanityCheck() {
            if (empty()) return false;
            auto repr = _contractions.at(0).u;
            return std::all_of(_contractions.begin(),_contractions.end(),[repr](Contraction c){return c.u == repr;});
        };

    public:
        ContractionIterator begin() {return _contractions.begin();};
        ContractionIterator end() {return _contractions.end();};
        auto cbegin() const {return _contractions.begin();};
        auto cend() const {return _contractions.end();};
        auto begin() const {return _contractions.begin();};
        auto end() const {return _contractions.end();};

        bool empty() const {return _contractions.empty();};
        uint64_t size() const {return _contractions.size();};

        Contraction at(int i) const {return _contractions.at(i);};

        ContractionGroup(std::initializer_list<Contraction> init) : _contractions(init) {
            ASSERT(sanityCheck());
            _representative = _contractions.front().u;
        };

        explicit ContractionGroup(std::vector<Contraction> init) : _contractions(std::move(init)) {
            ASSERT(sanityCheck());
            _representative = _contractions.front().u;
        }

        ContractionGroup(ContractionGroup& other) = default;
        ContractionGroup(const ContractionGroup& other) = default;

        HypernodeID getRepresentative() const {
            return _representative;
        }

        /// This function is linear in the number of contractions.
        bool contains(Contraction contraction) const {
            if (std::any_of(_contractions.begin(),_contractions.end(),[contraction](Contraction e){return e.u == contraction.u && e.v == contraction.v;})) {
                return true;
            }
            return false;
        }

        /// This function is quadratic in the group size.
        bool operator==(const ContractionGroup &rhs) const {
            bool rhsContainsThis = (std::all_of(_contractions.begin(),_contractions.end(),[rhs](Contraction e){return rhs.contains(e);}));
            bool thisContainsRhs = (std::all_of(rhs.begin(),rhs.end(),[this](Contraction e){return this->contains(e);}));
            return rhsContainsThis && thisContainsRhs;
        }

        /// This function is quadratic in the group size.
        bool operator!=(const ContractionGroup &rhs) const {
            return !(rhs == *this);
        }

        void debugPrint() const {
            std::cout << "\tGroup has " << size() << " contractions: \n";
            for (auto &c : *this) {
                std::cout << "\t\t(u: " << c.u << ", v: " << c.v << ")\n";
            }
        }
    };

    /// Adapts an ContractionIterator to an iterator over the contracted nodes in the contractions.
    /// I.e. this iterates over the 'v' value in the contractions that the given iterator defines.
    class ContractionToNodeIteratorAdaptor
        : public boost::iterator_adaptor<
                ContractionToNodeIteratorAdaptor,
                ContractionIterator,
                const HypernodeID>
        {

        using iterator = ContractionToNodeIteratorAdaptor;

        public:
            ContractionToNodeIteratorAdaptor() : ContractionToNodeIteratorAdaptor::iterator_adaptor_ () {}
            explicit ContractionToNodeIteratorAdaptor(const iterator::iterator_adaptor_::base_type& t) : ContractionToNodeIteratorAdaptor::iterator_adaptor_(t) {}

        private:
            friend class boost::iterator_core_access;

            // ! Adapting happens here: Dereference by redirecting to the contracted vertex
            const HypernodeID& dereference() const {return this->base()->v;}
    };

}

#endif //KAHYPAR_ASYNCH_COMMON_H

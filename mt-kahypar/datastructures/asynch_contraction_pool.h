//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_queue.h>
#include "gmock/gmock.h"

#include <utility>
#include "hypergraph_common.h"
#include <list>

namespace mt_kahypar::ds
{

    using Contraction = Memento;
    typedef std::vector<Contraction>::iterator ContractionIterator;

    class ContractionGroup {


private:

    std::vector<Contraction> contractions;

    /// Check if all contractions have the same representative
    bool sanityCheck() {
        auto repr = contractions.at(0).u;
        return std::all_of(contractions.begin(),contractions.end(),[repr](Contraction c){return c.u == repr;});
    };

public:

        bool empty() const {return contractions.empty();};
        uint64_t size() const {return contractions.size();};

        Contraction at(int i) const {return contractions.at(i);};

        ContractionGroup() = default;
        ContractionGroup(std::initializer_list<Contraction> init) : contractions(init) {
            ASSERT(sanityCheck());
        };

        explicit ContractionGroup(std::vector<Contraction> init) : contractions(std::move(init)) {
            ASSERT(sanityCheck());
        }

    /// This function is linear in the number of contractions. Only use for debugging!
    bool contains(Contraction contraction) const;

    /// This function is quadratic in the group size. Only use for debugging!
    bool operator==(const ContractionGroup &rhs) const;

    /// This function is quadratic in the group size. Only use for debugging!
    bool operator!=(const ContractionGroup &rhs) const;

    ContractionIterator begin() {return contractions.begin();};
    ContractionIterator end() {return contractions.end();};
    auto cbegin() const {return contractions.begin();};
    auto cend() const {return contractions.end();};
    auto begin() const {return contractions.begin();};
    auto end() const {return contractions.end();};

        void debugPrint() const;
    };

//    ContractionIterator begin(ContractionGroup& group){return group.begin();};
//    ContractionIterator end(ContractionGroup& group){return group.end();};
//    auto cbegin(const ContractionGroup& group) {return group.begin();};
//    auto cend(const ContractionGroup& group){return group.end();};
//    auto begin(const ContractionGroup& group){return group.begin();};
//    auto end(const ContractionGroup& group){return group.end();};

/// Pure virtual interface for a contraction pool
class IContractionPool {

public:
        virtual ~IContractionPool() = default;

        virtual uint64_t unsafe_size() const = 0;
        virtual bool empty() const = 0;

        /// Contains functions only for testing
        virtual bool contains(Contraction contraction) const = 0;
        virtual bool contains(const ContractionGroup& group) const = 0;

        /*!
         * Picks any contraction group from the pool. This explicitly does not define the way a group is chosen and is therefore not the same as picking a random group.
         * @return a contraction group.
         */
        virtual ContractionGroup pickAnyGroup() = 0;

//        /*!
//         * Picks a random contraction group from the pool.
//         * @return a contraction group.
//         */
//        ContractionGroup pickRandomGroup();

        /*!
         * Inserts a contraction as a single-contraction group implying it can be worked on independently from any other contraction.
         * @param contraction the contraction to be inserted as the only contraction in a group
         */
        virtual void insertContraction(Contraction contraction) = 0;

        /*!
         * Inserts a contraction group implying the contractions in this group have to be worked on collectively. All contractions in the group thus must have the same representative.
         * @param group the contraction group to insert.
         */
        virtual void insertContractionGroup(const ContractionGroup& group) = 0;
};

class SequentialContractionPool : public IContractionPool {

    private:
       std::list<ContractionGroup> groups;

    public:

        ~SequentialContractionPool() override = default;

        uint64_t unsafe_size() const override;
        bool empty() const override;

        /// Contains functions only for testing
        bool contains(Contraction contraction) const override;
        bool contains(const ContractionGroup& group) const override;

        /*!
         * Picks any contraction group from the pool. This explicitly does not define the way a group is chosen and is therefore not the same as picking a random group.
         * @return a contraction group.
         */
        ContractionGroup pickAnyGroup() override;

        /*!
         * Inserts a contraction as a single-contraction group implying it can be worked on independently from any other contraction.
         * @param contraction the contraction to be inserted as the only contraction in a group
         */
        void insertContraction(Contraction contraction) override;

        /*!
         * Inserts a contraction group implying the contractions in this group have to be worked on collectively. All contractions in the group thus must have the same representative.
         * @param group the contraction group to insert.
         */
        void insertContractionGroup(const ContractionGroup& group) override;

        void debugPrint() {
            std::cout << "Pool currently has " << unsafe_size() << " groups: \n";
            for (auto it = groups.begin(); it != groups.end(); ++it) {
                auto group = *it;
                group.debugPrint();
            }
        };

    };


class MockContractionPool : public IContractionPool {
    public:
        MOCK_METHOD(void,insertContraction, (Contraction),(override));
        MOCK_METHOD(void,insertContractionGroup, (const ContractionGroup&),(override));
        MOCK_METHOD(uint64_t, unsafe_size,(), (const,override));
        MOCK_METHOD(bool, empty,(), (const,override));
        MOCK_METHOD(bool, contains, (Contraction), (const,override));
        MOCK_METHOD(bool, contains, (const ContractionGroup&), (const,override));
        MOCK_METHOD(ContractionGroup, pickAnyGroup,(), (override));
};

} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H

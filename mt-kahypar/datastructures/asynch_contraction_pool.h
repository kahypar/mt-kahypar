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
#include "uncontraction_group_tree.h"
#include <list>

namespace mt_kahypar::ds
{


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

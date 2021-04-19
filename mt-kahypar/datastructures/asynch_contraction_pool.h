//
// Created by mlaupichler on 19.04.21.
//

#ifndef KAHYPAR_ASYNCH_CONTRACTION_POOL_H
#define KAHYPAR_ASYNCH_CONTRACTION_POOL_H

#include <tbb/concurrent_queue.h>
#include "hypergraph_common.h"

namespace mt_kahypar::ds {

    using Contraction = Memento;
    typedef std::vector<Contraction>::iterator ContractionIterator;

    class ContractionGroup {


private:

    std::vector<Contraction> contractions;

public:

        bool empty() const {return contractions.empty();};
        uint64_t size() const {return contractions.size();};

        Contraction at(int i) const {return contractions.at(i);};

        ContractionGroup() = default;
        ContractionGroup(std::initializer_list<Contraction> init) : contractions(init) {
            //assert that they all have the same representative
            auto repr = contractions.at(0).u;
            for(auto &contr : contractions) {
                ASSERT(contr.u == repr);
            }
        };

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

class AsynchContractionPool {

    private:
        tbb::concurrent_queue<ContractionGroup> queue;

    public:

        uint64_t unsafe_size();
        bool empty();

        /// Contains functions only for testing
        bool contains(Contraction contraction);
        bool contains(const ContractionGroup& group);

        /*!
         * Picks any contraction group from the pool. This explicitly does not define the way a group is chosen and is therefore not the same as picking a random group.
         * @return a contraction group.
         */
        ContractionGroup pickAnyGroup();

//        /*!
//         * Picks a random contraction group from the pool.
//         * @return a contraction group.
//         */
//        ContractionGroup pickRandomGroup();

        /*!
         * Inserts a contraction as a single-contraction group implying it can be worked on independently from any other contraction.
         * @param contraction the contraction to be inserted as the only contraction in a group
         */
        void insertContraction(Contraction contraction);

        /*!
         * Inserts a contraction group implying the contractions in this group have to be worked on collectively. All contractions in the group thus must have the same representative.
         * @param group the contraction group to insert.
         */
        void insertContractionGroup(const ContractionGroup& group);

        void debugPrintPool() {
            std::cout << "Pool currently has " << unsafe_size() << " groups: \n";
            for (auto it = queue.unsafe_begin(); it != queue.unsafe_end(); ++it) {
                ContractionGroup group = *it;
                group.debugPrint();
            }
        };

    };

} // namespace mt_kahypar
#endif //KAHYPAR_ASYNCH_CONTRACTION_POOL_H

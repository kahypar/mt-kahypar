//
// Created by mlaupichler on 19.04.21.
//

#include "asynch_contraction_pool.h"

uint64_t mt_kahypar::ds::SequentialContractionPool::unsafe_size() const {
    return groups.size();
}

bool mt_kahypar::ds::SequentialContractionPool::empty() const {
    return groups.empty();
}

void mt_kahypar::ds::SequentialContractionPool::insertContractionGroup(const mt_kahypar::ds::ContractionGroup& group) {
    ASSERT(!group.empty());

    groups.push_back(group);
}

void mt_kahypar::ds::SequentialContractionPool::insertContraction(mt_kahypar::ds::Contraction contraction) {
    ContractionGroup group = {contraction};
    insertContractionGroup(group);
}

bool mt_kahypar::ds::SequentialContractionPool::contains(const mt_kahypar::ds::ContractionGroup& group) const {

    std::cout << "Searching for: \n";
    group.debugPrint();

    for (auto it = groups.begin(); it != groups.end(); ++it) {
//        std::cout << it->size() << "\n";
        ContractionGroup poolEl = *it;
        if (poolEl == group) {
            std::cout << "Found in pool: \n";
            poolEl.debugPrint();
            std::cout << "Searched for: \n";
            group.debugPrint();
            return true;
        }
    }
    return false;
}

bool mt_kahypar::ds::SequentialContractionPool::contains(mt_kahypar::ds::Contraction contraction) const {
    for (auto it = groups.begin(); it != groups.end(); ++it) {
        ContractionGroup group = *it;
        if (group.contains(contraction)) return true;
    }
    return false;
}


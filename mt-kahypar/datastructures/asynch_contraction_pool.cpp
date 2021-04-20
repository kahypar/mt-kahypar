//
// Created by mlaupichler on 19.04.21.
//

#include "asynch_contraction_pool.h"

uint64_t mt_kahypar::ds::AsynchContractionPool::unsafe_size() {
    return groups.unsafe_size();
}

bool mt_kahypar::ds::AsynchContractionPool::empty() {
    return groups.empty();
}

void mt_kahypar::ds::AsynchContractionPool::insertContractionGroup(const mt_kahypar::ds::ContractionGroup& group) {
    ASSERT(!group.empty());

    groups.push(group);
}

void mt_kahypar::ds::AsynchContractionPool::insertContraction(mt_kahypar::ds::Contraction contraction) {
    ContractionGroup group = {contraction};
    insertContractionGroup(group);
}

bool mt_kahypar::ds::AsynchContractionPool::contains(const mt_kahypar::ds::ContractionGroup& group) {

    std::cout << "Searching for: \n";
    group.debugPrint();

    for (auto it = groups.unsafe_begin(); it != groups.unsafe_end(); ++it) {
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

bool mt_kahypar::ds::AsynchContractionPool::contains(mt_kahypar::ds::Contraction contraction) {
    for (auto it = groups.unsafe_begin(); it != groups.unsafe_end(); ++it) {
        ContractionGroup group = *it;
        if (group.contains(contraction)) return true;
    }
    return false;
}

mt_kahypar::ds::ContractionGroup mt_kahypar::ds::AsynchContractionPool::pickAnyGroup() {
    ASSERT(!groups.empty());
    ContractionGroup result;
    groups.try_pop(result);
    return result;
}


bool mt_kahypar::ds::ContractionGroup::contains(mt_kahypar::ds::Contraction contraction) const {
    if (std::any_of(contractions.begin(),contractions.end(),[contraction](Contraction e){return e.u == contraction.u && e.v == contraction.v;})) {
        return true;
    }
    return false;
}

bool mt_kahypar::ds::ContractionGroup::operator==(const mt_kahypar::ds::ContractionGroup &rhs) const {

    bool rhsContainsThis = (std::all_of(contractions.begin(),contractions.end(),[rhs](Contraction e){return rhs.contains(e);}));
    bool thisContainsRhs = (std::all_of(rhs.begin(),rhs.end(),[this](Contraction e){return this->contains(e);}));

    return rhsContainsThis && thisContainsRhs;
}

bool mt_kahypar::ds::ContractionGroup::operator!=(const mt_kahypar::ds::ContractionGroup &rhs) const {
    return !(rhs == *this);
}

void mt_kahypar::ds::ContractionGroup::debugPrint() const {
    std::cout << "\tGroup has " << size() << " contractions: \n";
    for (auto &c : *this) {
        std::cout << "\t\t(u: " << c.u << ", v: " << c.v << ")\n";
    }
}

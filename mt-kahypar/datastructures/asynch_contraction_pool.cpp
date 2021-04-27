//
// Created by mlaupichler on 19.04.21.
//

#include "asynch_contraction_pool.h"

uint32_t mt_kahypar::ds::SequentialContractionGroupPool::getNumActive() const {
    return _active.size();
}

const mt_kahypar::ds::ContractionGroup &
mt_kahypar::ds::SequentialContractionGroupPool::group(mt_kahypar::ds::ContractionGroupID id) const {
    return _hierarchy->group(id);
}

mt_kahypar::ds::ContractionGroupID mt_kahypar::ds::SequentialContractionGroupPool::pickAnyActiveID() {
    ASSERT(!_active.empty());
    // This implementation is LIFO
    auto picked = _active.back();
    _active.pop_back();
    return picked;
}

void mt_kahypar::ds::SequentialContractionGroupPool::activateSuccessors(mt_kahypar::ds::ContractionGroupID id) {

    auto succs = _hierarchy->successors(id);
    for (auto s : succs) {
        _active.push_back(s);
    }
}

bool mt_kahypar::ds::SequentialContractionGroupPool::hasActive() const {
    return !(_active.empty());
}

mt_kahypar::ds::BlockedGroupIDIterator mt_kahypar::ds::SequentialContractionGroupPool::all() const {
    return _hierarchy->all();
}

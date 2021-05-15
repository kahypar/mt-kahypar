//
// Created by mlaupichler on 15.05.21.
//

#ifndef KAHYPAR_MOCK_GROUP_HIERARCHY_H
#define KAHYPAR_MOCK_GROUP_HIERARCHY_H

#include "mt-kahypar/datastructures/asynch/asynch_common.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace mt_kahypar::ds{

    class MockGroupHierarchy {
    public:
        MOCK_METHOD(size_t,getVersion,(), (const));
        MOCK_METHOD(uint32_t, getNumGroups,(),(const));
        MOCK_METHOD(const ContractionGroup&, group,(ContractionGroupID), (const));
        MOCK_METHOD(ContractionGroupID, predecessor,(ContractionGroupID), (const));
        MOCK_METHOD(ContractionGroupIDIteratorRange,successors,(ContractionGroupID),(const));
        MOCK_METHOD(ContractionGroupIDIteratorRange,roots,(),(const));
        MOCK_METHOD(BlockedGroupIDIterator,all,(),(const));
    };

}




#endif //KAHYPAR_MOCK_GROUP_HIERARCHY_H

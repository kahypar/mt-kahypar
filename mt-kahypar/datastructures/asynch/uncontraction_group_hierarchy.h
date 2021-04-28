//
// Created by mlaupichler on 28.04.21.
//

#ifndef KAHYPAR_UNCONTRACTION_GROUP_HIERARCHY_H
#define KAHYPAR_UNCONTRACTION_GROUP_HIERARCHY_H

#include "gmock/gmock.h"

#include "asynch_common.h"

namespace mt_kahypar::ds {

    /// Pure virtual interface to define functions for a hierarchy of UncontractionGroups in a version.
    /// In this hierarchy an UncontractionGroup is supposed to be uncontracted only once its predecessor has been uncontracted.
    /// This expresses the minimal order between uncontraction groups based on the order of contractions.
    class IUncontractionGroupHierarchy {

    public:

        static constexpr size_t kInvalidVersion = std::numeric_limits<size_t>::max();

        virtual size_t getVersion() const = 0;
        virtual uint32_t getNumGroups() const = 0;
        virtual const ContractionGroup& group(ContractionGroupID id) const = 0;
        virtual ContractionGroupID predecessor(ContractionGroupID id) const = 0;
        virtual ContractionGroupIDIteratorRange successors(ContractionGroupID id) const = 0;
        virtual ContractionGroupIDIteratorRange roots() const = 0;

        virtual BlockedGroupIDIterator all() const = 0;

        virtual ~IUncontractionGroupHierarchy() = default;
    };

    class MockGroupHierarchy : public IUncontractionGroupHierarchy {
    public:
        MOCK_METHOD(size_t,getVersion,(), (const,override));
        MOCK_METHOD(uint32_t, getNumGroups,(),(const, override));
        MOCK_METHOD(const ContractionGroup&, group,(ContractionGroupID), (const, override));
        MOCK_METHOD(ContractionGroupID, predecessor,(ContractionGroupID), (const, override));
        MOCK_METHOD(ContractionGroupIDIteratorRange,successors,(ContractionGroupID),(const, override));
        MOCK_METHOD(ContractionGroupIDIteratorRange,roots,(),(const, override));
        MOCK_METHOD(BlockedGroupIDIterator,all,(),(const, override));
    };



}

#endif //KAHYPAR_UNCONTRACTION_GROUP_HIERARCHY_H

//
// Created by mlaupichler on 24.04.21.
//

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/asynch/uncontraction_group_tree.h"

namespace mt_kahypar::ds {

    TEST(AUncontractionGroupTree,SuccessorsGivesEmptyRangeWhenNoSuccessors) {
        size_t version = 0;

        ContractionTree contractionTree;
        contractionTree.initialize(2);
        contractionTree.setInterval(0, 0, 1);
        contractionTree.setParent(1, 0, version); contractionTree.setInterval(1, 2, 3);

        contractionTree.finalize(1);

        ContractionGroup expectedRootGroup1 = { Contraction {0, 1}};

        auto contractionTreeCopy = contractionTree.copy();
        UncontractionGroupTree groupTree = UncontractionGroupTree(contractionTreeCopy, version);

        ASSERT(groupTree.getNumGroups() == 1);
        ASSERT(groupTree.getVersion() == version);

        auto roots = groupTree.roots();
        ASSERT(std::distance(roots.begin(),roots.end()) == 1);
        auto rootID1 = *roots.begin();
        const auto& rootGroup1 = groupTree.group(rootID1);
        ASSERT(rootGroup1 == expectedRootGroup1);

        auto succs = groupTree.successors(rootID1);
        ASSERT(succs.empty());
        for (auto s : succs) {
            ASSERT(false, "Can iterate in empty successors range!");
        }

    }

    TEST(AUncontractionGroupTree,CreatesCorrectTreeForUniformVersion1) {

        size_t version = 0;

        ContractionTree contractionTree;
        contractionTree.initialize(4);
//    contractionTree.setParent(0, 0, version);
        contractionTree.setInterval(0, 0, 1);
        contractionTree.setParent(1, 0, version); contractionTree.setInterval(1, 2, 3); // 1,2,3 are siblings but only 2 and 3 have overlapping intervals
        contractionTree.setParent(2, 0, version); contractionTree.setInterval(2, 4, 5);
        contractionTree.setParent(3, 0, version); contractionTree.setInterval(3, 5, 6);
//    contractionTree.setParent(4, 1, version); contractionTree.setInterval(4,5,6); // 4,5 are siblings but do not have overlapping intervals
//    contractionTree.setParent(5, 1, version); contractionTree.setInterval(5,7,8);
//    contractionTree.setParent(6, 3, version); contractionTree.setInterval(6,9,10);

        contractionTree.finalize(1);

        ContractionGroup expectedRootGroup1 = { Contraction {0, 2}, Contraction {0, 3}};
        ContractionGroup expectedSuccGroup1 = { Contraction {0, 1}};

        auto contractionTreeCopy = contractionTree.copy();
        UncontractionGroupTree groupTree = UncontractionGroupTree(contractionTreeCopy, version);

        ASSERT(groupTree.getNumGroups() == 2);
        ASSERT(groupTree.getVersion() == version);

        auto roots = groupTree.roots();
        ASSERT(std::distance(roots.begin(),roots.end()) == 1);
        auto rootID1 = *roots.begin();
        const auto& rootGroup1 = groupTree.group(rootID1);
        ASSERT(rootGroup1 == expectedRootGroup1);

        auto succsOfRoot1 = groupTree.successors(rootID1);
        ASSERT(std::distance(succsOfRoot1.begin(),succsOfRoot1.end()) == 1);
        auto succID1 = *succsOfRoot1.begin();
        const auto& succGroup1 = groupTree.group(succID1);
        ASSERT(succGroup1 == expectedSuccGroup1);
        ASSERT(groupTree.predecessor(succID1) == rootID1);

    }

    TEST(AUncontractionGroupTree,CreatesCorrectTreeForUniformVersion2) {

        size_t version = 0;

        // roots 0 and 4
        ContractionTree contractionTree;
        contractionTree.initialize(7);
//    contractionTree.setParent(0, 0, version);
        contractionTree.setInterval(0, 0, 1);
        contractionTree.setParent(1, 0, version); contractionTree.setInterval(1, 2, 3); // 1,2,3 are siblings but only 2 and 3 have overlapping intervals
        contractionTree.setParent(2, 0, version); contractionTree.setInterval(2, 4, 5);
        contractionTree.setParent(3, 0, version); contractionTree.setInterval(3, 5, 6);
//    contractionTree.setParent(4, 1, version);
        contractionTree.setInterval(4, 5, 6);
        contractionTree.setParent(5, 4, version); contractionTree.setInterval(5, 7, 8); // Root group does not have horizontal children but one vertical child
        contractionTree.setParent(6, 5, version); contractionTree.setInterval(6, 9, 10);

        contractionTree.finalize(1);

        ContractionGroup expectedRootGroup1 = { Contraction {0, 2}, Contraction {0, 3}};
        ContractionGroup expectedSuccGroup1 = { Contraction {0, 1}};
        ContractionGroup expectedRootGroup2 = {Contraction {4, 5}};
        ContractionGroup expectedSuccGroup2 = {Contraction {5, 6}};

        auto contractionTreeCopy = contractionTree.copy();
        UncontractionGroupTree groupTree = UncontractionGroupTree(contractionTreeCopy, version);

        ASSERT(groupTree.getNumGroups() == 4);
        ASSERT(groupTree.getVersion() == version);

        auto roots = groupTree.roots();
        ASSERT(std::distance(roots.begin(),roots.end()) == 2);

        bool seenExpected1 = false;
        bool seenExpected2 = false;

        for (auto r : roots) {
            auto rootID = r;
            const auto& rootGroup = groupTree.group(rootID);
            ASSERT(rootGroup == expectedRootGroup1 || rootGroup == expectedRootGroup2);

            if (rootGroup == expectedRootGroup1) {
                seenExpected1 = true;
                auto succsOfRoot1 = groupTree.successors(rootID);
                ASSERT(std::distance(succsOfRoot1.begin(),succsOfRoot1.end()) == 1);
                auto succID1 = *succsOfRoot1.begin();
                const auto& succGroup1 = groupTree.group(succID1);
                ASSERT(succGroup1 == expectedSuccGroup1);
                ASSERT(groupTree.predecessor(succID1) == rootID);
            }
            if (rootGroup == expectedRootGroup2) {
                seenExpected2 = true;
                auto succsOfRoot2 = groupTree.successors(rootID);
                ASSERT(std::distance(succsOfRoot2.begin(),succsOfRoot2.end()) == 1);
                auto succID2 = *succsOfRoot2.begin();
                const auto& succGroup2 = groupTree.group(succID2);
                ASSERT(succGroup2 == expectedSuccGroup2);
                ASSERT(groupTree.predecessor(succID2) == rootID);
            }
        }
        ASSERT(seenExpected1 && seenExpected2);
    }

    TEST(AUncontractionGroupTree,CreatesCorrectTreeForDifferentVersions1) {

        size_t version0 = 1;
        size_t version1 = 0;

        // roots 0 and 4
        ContractionTree contractionTree;
        contractionTree.initialize(7);
//    contractionTree.setParent(0, 0, version);
        contractionTree.setInterval(0, 0, 1);
        contractionTree.setParent(1, 0, version0); contractionTree.setInterval(1, 2, 3); // 1,2,3 are siblings but only 2 and 3 have overlapping intervals but 2 and 3 have differing versions!
        contractionTree.setParent(2, 0, version0); contractionTree.setInterval(2, 4, 5);
        contractionTree.setParent(3, 0, version1); contractionTree.setInterval(3, 4, 6);
        contractionTree.setParent(4, 1, version1); contractionTree.setInterval(4, 0, 1); // 4,5 are siblings in the same version and overlap
        contractionTree.setParent(5, 1, version1); contractionTree.setInterval(5, 0, 1);
        contractionTree.setParent(6, 3, version1); contractionTree.setInterval(6, 1, 2);

        contractionTree.finalize(2);

        // Version 0 test
        ContractionGroup expectedVersion0RootGroup1 = { Contraction {0, 2}};
        ContractionGroup expectedVersion0SuccGroup1 = { Contraction {0, 1}};

        auto contractionTreeCopy = contractionTree.copy();
        UncontractionGroupTree groupTreeVersion0 = UncontractionGroupTree(contractionTreeCopy, version0);

        ASSERT(groupTreeVersion0.getNumGroups() == 2);
        ASSERT(groupTreeVersion0.getVersion() == version0);

        auto rootsVersion0 = groupTreeVersion0.roots();
        ASSERT(std::distance(rootsVersion0.begin(),rootsVersion0.end()) == 1);

        auto rootIDVersion0 = *rootsVersion0.begin();
        const auto& rootGroupVersion0 = groupTreeVersion0.group(rootIDVersion0);
        ASSERT(rootGroupVersion0 == expectedVersion0RootGroup1);
        auto succsOfRoot1Version0 = groupTreeVersion0.successors(rootIDVersion0);
        ASSERT(std::distance(succsOfRoot1Version0.begin(),succsOfRoot1Version0.end()) == 1);
        auto succID1 = *succsOfRoot1Version0.begin();
        const auto& succGroup1 = groupTreeVersion0.group(succID1);
        ASSERT(succGroup1 == expectedVersion0SuccGroup1);
        ASSERT(groupTreeVersion0.predecessor(succID1) == rootIDVersion0);

        // Version 1 test
        ContractionGroup expectedVersion1RootGroup1 = { Contraction {0, 3}};
        ContractionGroup expectedVersion1RootGroup2 = { Contraction {1, 4}, {1, 5}};
        ContractionGroup expectedVersion1SuccGroup1 = { Contraction {3, 6}};

        auto contractionTreeCopy2 = contractionTree.copy();
        UncontractionGroupTree groupTreeVersion1 = UncontractionGroupTree(contractionTreeCopy2, version1);

        ASSERT(groupTreeVersion1.getNumGroups() == 3);
        ASSERT(groupTreeVersion1.getVersion() == version1);

        auto roots = groupTreeVersion1.roots();
        ASSERT(std::distance(roots.begin(),roots.end()) == 2);

        bool seenExpected1 = false;
        bool seenExpected2 = false;

        for (auto r : roots) {
            auto rootID = r;
            const auto& rootGroup = groupTreeVersion1.group(rootID);
            ASSERT(rootGroup == expectedVersion1RootGroup1 || rootGroup == expectedVersion1RootGroup2);

            if (rootGroup == expectedVersion1RootGroup1) {
                seenExpected1 = true;
                auto succsOfRoot1 = groupTreeVersion1.successors(rootID);
                ASSERT(std::distance(succsOfRoot1.begin(),succsOfRoot1.end()) == 1);
                auto succID1 = *succsOfRoot1.begin();
                const auto& succGroup1 = groupTreeVersion1.group(succID1);
                ASSERT(succGroup1 == expectedVersion1SuccGroup1);
                ASSERT(groupTreeVersion1.predecessor(succID1) == rootID);
            }
            if (rootGroup == expectedVersion1RootGroup2) {
                seenExpected2 = true;
                auto succsOfRoot2 = groupTreeVersion1.successors(rootID);
                ASSERT(std::distance(succsOfRoot2.begin(),succsOfRoot2.end()) == 0);
            }
        }
        ASSERT(seenExpected1 && seenExpected2);

    }
} // namespace mt_kahypar::ds



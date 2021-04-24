//
// Created by mlaupichler on 24.04.21.
//

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/uncontraction_group_tree.h"

namespace mt_kahypar::ds {

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
        auto succGroup1 = groupTree.group(succID1);
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

        for (auto r : roots) {
            auto rootID = r;
            auto rootGroup = groupTree.group(rootID);
            ASSERT(rootGroup == expectedRootGroup1 || rootGroup == expectedRootGroup2);

            if (rootGroup == expectedRootGroup1) {
                auto succsOfRoot1 = groupTree.successors(rootID);
                ASSERT(std::distance(succsOfRoot1.begin(),succsOfRoot1.end()) == 1);
                auto succID1 = *succsOfRoot1.begin();
                auto succGroup1 = groupTree.group(succID1);
                ASSERT(succGroup1 == expectedSuccGroup1);
                ASSERT(groupTree.predecessor(succID1) == rootID);
            }
            if (rootGroup == expectedRootGroup2) {
                auto succsOfRoot2 = groupTree.successors(rootID);
                ASSERT(std::distance(succsOfRoot2.begin(),succsOfRoot2.end()) == 1);
                auto succID2 = *succsOfRoot2.begin();
                auto succGroup2 = groupTree.group(succID2);
                ASSERT(succGroup2 == expectedSuccGroup2);
                ASSERT(groupTree.predecessor(succID2) == rootID);
            }
        }
    }

    TEST(AContractionGroup, EqualityTest) {
        Contraction c1 = {1,2};
        Contraction c2 = {1,4};
        Contraction c3 = {1,3};
        ContractionGroup group1 = {c1,c2};
        ContractionGroup group2 = {c1,c3};
        ContractionGroup group3 = {c2,c1};

        ASSERT_FALSE(group1 == group2);
        ASSERT_TRUE(group1 == group3);

        ContractionGroup group4 = {c1};
        ASSERT_FALSE(group1 == group4);
        ASSERT_FALSE(group4 == group1);
    }

    TEST(AContractionGroupDeathTest, ConstructFailsIfEmptyDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";

        EXPECT_DEATH(ContractionGroup group = {},"");
        std::vector<Contraction> emptyVec;
        EXPECT_DEATH(ContractionGroup group(emptyVec),"");
    }

    TEST(AContractionGroupDeathTest, ConstructionFailsIfGroupHasDifferentRepsDeathTest) {
        testing::FLAGS_gtest_death_test_style="threadsafe";
        auto pool = SequentialContractionPool();
        Contraction c1 = {1,2};
        Contraction c2 = {2,3};
        ASSERT_DEATH(ContractionGroup({c1,c2}),"");
    }
} // namespace mt_kahypar::ds



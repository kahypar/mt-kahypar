//
// Created by mlaupichler on 25.07.21.
//

#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"

#include "mt-kahypar/datastructures/async/node_region_comparator.h"

namespace mt_kahypar::ds {

    using Fac = DynamicHypergraphFactory;
    using HG = DynamicHypergraph;

    TEST(ANodeRegionComparator, SimilarityAlwaysZeroForDegreeZeroNodes) {

      HG hg = Fac::construct(3, 1, {{2}});

      auto comp = NodeRegionComparator<HG>(1);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 1), 0.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 2), 0.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 2), 0.0);
    }

    TEST(ANodeRegionComparator, SimilarityAlwaysOnForSameNode) {

      HG hg = Fac::construct(3, 1, {{1,2}});

      auto comp = NodeRegionComparator<HG>(1);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 0), 1.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 1), 1.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(2, 2), 1.0);
    }

    TEST(ANodeRegionComparator, SimilarityAlwaysOneForEquivalentNodes) {

      HG hg = Fac::construct(6, 3, {{0,1,2,3,4}, {3,4}, {3,4,5}});

      auto comp = NodeRegionComparator<HG>(3);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 1), 1.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 2), 1.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 2), 1.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(3, 4), 1.0);
    }

    TEST(ANodeRegionComparator, SimilarityAlwaysZeroForNodesWithoutIntersection) {

      HG hg = Fac::construct(4, 2, {{0,1}, {2,3}});

      auto comp = NodeRegionComparator<HG>(2);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 2), 0.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 3), 0.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 2), 0.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 3), 0.0);
    }

    TEST(ANodeRegionComparator, SimilarityWhenOneSignatureIsSubsumed) {

      HG hg = Fac::construct(2, 10,
                             {{0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {1}, {1}, {1}, {1}, {1}});

      auto comp = NodeRegionComparator<HG>(10);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 1), 0.5);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 0), 0.5);
    }

    TEST(ANodeRegionComparator, SimilarityBasedOnlyOnLargestEdges) {

      // A larger edge for both 0 and 1 will make up their entire signatures and they should have similarity zero even
      // though they are actually neighbors in a smaller edge
      HG hg = Fac::construct(4, 3,
                             {{0,1}, {0,2,3}, {1,2,3}});

      auto comp = NodeRegionComparator<HG>(1);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 1), 0.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 0), 0.0);
    }

    TEST(ANodeRegionComparator, GeneralSimilarity1) {

      HG hg = Fac::construct(4, 3,
                             {{0,1}, {0,2,3}, {1,2,3}});

      auto comp = NodeRegionComparator<HG>(3);
      comp.calculateSignaturesParallel(&hg);

      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 1), ((double) 1.0 / 3.0));
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 2), ((double) 1.0 / 3.0));
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(0, 3), ((double) 1.0 / 3.0));
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 2), ((double) 1.0 / 3.0));
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(1, 3), ((double) 1.0 / 3.0));
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(2, 3), 1.0);
      ASSERT_DOUBLE_EQ(comp.regionSimilarity(2, 3), 1.0);
    }

    TEST(ANodeRegionComparator, GeneralSimilarity2) {

      HG hg = Fac::construct(2, 10,
                             {{0,1}, {0}, {0}, {0}, {0}, {1}, {1}, {1}, {1}, {1}});

      auto comp1 = NodeRegionComparator<HG>(2);
      comp1.calculateSignaturesParallel(&hg);
      ASSERT_DOUBLE_EQ(comp1.regionSimilarity(0, 1), ((double) 1.0 / 3.0));

      auto comp2 = NodeRegionComparator<HG>(10);
      comp2.calculateSignaturesParallel(&hg);
      ASSERT_DOUBLE_EQ(comp2.regionSimilarity(0, 1), 0.1);
    }

} // end namespace mt_kahypar::ds

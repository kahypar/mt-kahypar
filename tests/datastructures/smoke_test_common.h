//
// Created by mlaupichler on 30.04.21.
//

#ifndef KAHYPAR_SMOKE_TEST_COMMON_H
#define KAHYPAR_SMOKE_TEST_COMMON_H

#include "mt-kahypar/datastructures/dynamic_hypergraph.h"

namespace mt_kahypar {
    namespace ds {

        using DynamicPartitionedHypergraph = PartitionedHypergraph<DynamicHypergraph, DynamicHypergraphFactory, GainCache>;

        void verifyEqualityOfHypergraphs(const DynamicHypergraph& expected_hypergraph,
                                         const DynamicHypergraph& actual_hypergraph);

        HyperedgeWeight compute_km1(DynamicPartitionedHypergraph& partitioned_hypergraph);

        void verifyGainCache(DynamicPartitionedHypergraph& partitioned_hypergraph);

        void verifyNumIncidentCutHyperedges(const DynamicPartitionedHypergraph& partitioned_hypergraph);

        DynamicHypergraph generateRandomHypergraph(const HypernodeID num_hypernodes,
                                                   const HyperedgeID num_hyperedges,
                                                   const HypernodeID max_edge_size);

        BatchVector generateRandomContractions(const HypernodeID num_hypernodes,
                                               const HypernodeID num_contractions,
                                               const bool multi_versioned = true);

        void generateRandomPartition(DynamicPartitionedHypergraph& partitioned_hypergraph);

    } // namespace ds
} // namespace mt_kahypar

#endif //KAHYPAR_SMOKE_TEST_COMMON_H

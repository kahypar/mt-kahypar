//
// Created by mlaupichler on 30.04.21.
//

#include "smoke_test_common.h"
#include "gtest/gtest.h"

namespace mt_kahypar {
    namespace ds {

        void verifyEqualityOfHypergraphs(const DynamicHypergraph& expected_hypergraph,
                                         const DynamicHypergraph& actual_hypergraph) {
            parallel::scalable_vector<HyperedgeID> expected_incident_edges;
            parallel::scalable_vector<HyperedgeID> actual_incident_edges;
            for ( const HypernodeID& hn : expected_hypergraph.nodes() ) {
                ASSERT_TRUE(actual_hypergraph.nodeIsEnabled(hn));
                ASSERT_EQ(expected_hypergraph.nodeWeight(hn), actual_hypergraph.nodeWeight(hn));
                ASSERT_EQ(expected_hypergraph.nodeDegree(hn), actual_hypergraph.nodeDegree(hn));
                for ( const HyperedgeID he : expected_hypergraph.incidentEdges(hn) ) {
                    expected_incident_edges.push_back(he);
                }
                for ( const HyperedgeID he : actual_hypergraph.incidentEdges(hn) ) {
                    actual_incident_edges.push_back(he);
                }
                std::sort(expected_incident_edges.begin(), expected_incident_edges.end());
                std::sort(actual_incident_edges.begin(), actual_incident_edges.end());
                ASSERT_EQ(expected_incident_edges.size(), actual_incident_edges.size());
                for ( size_t i = 0; i < expected_incident_edges.size(); ++i ) {
                    ASSERT_EQ(expected_incident_edges[i], actual_incident_edges[i]);
                }
                expected_incident_edges.clear();
                actual_incident_edges.clear();
            }

            parallel::scalable_vector<HypernodeID> expected_pins;
            parallel::scalable_vector<HypernodeID> actual_pins;
            for ( const HyperedgeID& he : expected_hypergraph.edges() ) {
                for ( const HyperedgeID he : expected_hypergraph.pins(he) ) {
                    expected_pins.push_back(he);
                }
                for ( const HyperedgeID he : actual_hypergraph.pins(he) ) {
                    actual_pins.push_back(he);
                }
                std::sort(expected_pins.begin(), expected_pins.end());
                std::sort(actual_pins.begin(), actual_pins.end());
                ASSERT_EQ(expected_pins.size(), actual_pins.size());
                for ( size_t i = 0; i < expected_pins.size(); ++i ) {
                    ASSERT_EQ(expected_pins[i], actual_pins[i]);
                }
                expected_pins.clear();
                actual_pins.clear();
            }
        }

        HyperedgeWeight compute_km1(DynamicPartitionedHypergraph& partitioned_hypergraph) {
            HyperedgeWeight km1 = 0;
            for (const HyperedgeID& he : partitioned_hypergraph.edges()) {
                km1 += std::max(partitioned_hypergraph.connectivity(he) - 1, 0) * partitioned_hypergraph.edgeWeight(he);
            }
            return km1;
        }

        void verifyGainCache(DynamicPartitionedHypergraph& partitioned_hypergraph) {

            ASSERT_TRUE(partitioned_hypergraph.checkTrackedPartitionInformation());

            const PartitionID k = partitioned_hypergraph.k();
            utils::Randomize& rand = utils::Randomize::instance();
            HyperedgeWeight km1_before = compute_km1(partitioned_hypergraph);
            HyperedgeWeight expected_gain = 0;
            for ( const HypernodeID& hn : partitioned_hypergraph.nodes() ) {
                const PartitionID from = partitioned_hypergraph.partID(hn);
                PartitionID to = rand.getRandomInt(0, k - 1, sched_getcpu());
                if ( from == to ) to = (to + 1) % k;
                expected_gain += partitioned_hypergraph.km1Gain(hn, from, to);
                partitioned_hypergraph.changeNodePartWithGainCacheUpdate(hn, from, to);
            }
            HyperedgeWeight km1_after = compute_km1(partitioned_hypergraph);
            ASSERT_EQ(expected_gain, km1_before - km1_after) << V(expected_gain) << V(km1_before) << V(km1_after);

        }

        void verifyNumIncidentCutHyperedges(const DynamicPartitionedHypergraph& partitioned_hypergraph) {
            partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
                HypernodeID expected_num_cut_hyperedges = 0;
                for ( const HyperedgeID& he : partitioned_hypergraph.incidentEdges(hn) ) {
                    if ( partitioned_hypergraph.connectivity(he) > 1 ) {
                        ++expected_num_cut_hyperedges;
                    }
                }
                ASSERT_EQ(expected_num_cut_hyperedges, partitioned_hypergraph.numIncidentCutHyperedges(hn));
            });
        }

        DynamicHypergraph generateRandomHypergraph(const HypernodeID num_hypernodes,
                                                   const HyperedgeID num_hyperedges,
                                                   const HypernodeID max_edge_size) {
            parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> hyperedges;
            utils::Randomize& rand = utils::Randomize::instance();
            for ( size_t i = 0; i < num_hyperedges; ++i ) {
                parallel::scalable_vector<HypernodeID> net;
                const size_t edge_size = rand.getRandomInt(2, max_edge_size, sched_getcpu());
                for ( size_t i = 0; i < edge_size; ++i ) {
                    const HypernodeID pin = rand.getRandomInt(0, num_hypernodes - 1, sched_getcpu());
                    if ( std::find(net.begin(), net.end(), pin) == net.end() ) {
                        net.push_back(pin);
                    }
                }
                hyperedges.emplace_back(std::move(net));
            }
            return DynamicHypergraphFactory::construct(num_hypernodes, num_hyperedges, hyperedges);
        }

        BatchVector generateRandomContractions(const HypernodeID num_hypernodes,
                                               const HypernodeID num_contractions,
                                               const bool multi_versioned) {
            ASSERT(num_contractions < num_hypernodes);
            HypernodeID tmp_num_contractions = num_contractions;
            BatchVector contractions;
            parallel::scalable_vector<HypernodeID> active_hns(num_hypernodes);
            std::iota(active_hns.begin(), active_hns.end(), 0);
            utils::Randomize& rand = utils::Randomize::instance();
            const int cpu_id = sched_getcpu();
            while ( tmp_num_contractions > 0 ) {
                HypernodeID current_num_contractions = tmp_num_contractions;
                if ( multi_versioned && current_num_contractions > 25 ) current_num_contractions /= 2;
                contractions.emplace_back();
                for ( size_t i = 0; i < current_num_contractions; ++i ) {
                    ASSERT(active_hns.size() >= 2UL);
                    int idx_1 = rand.getRandomInt(0, static_cast<int>(active_hns.size() - 1), cpu_id);
                    int idx_2 = rand.getRandomInt(0, static_cast<int>(active_hns.size() - 1), cpu_id);
                    if ( idx_1 == idx_2 ) {
                        idx_2 = (idx_2 + 1) % active_hns.size();
                    }
                    contractions.back().push_back(Memento { active_hns[idx_1], active_hns[idx_2] });
                    std::swap(active_hns[idx_2], active_hns.back());
                    active_hns.pop_back();
                }
                tmp_num_contractions -= current_num_contractions;
            }
            return contractions;
        }

        void generateRandomPartition(DynamicPartitionedHypergraph& partitioned_hypergraph) {
            const PartitionID k = partitioned_hypergraph.k();
            utils::Randomize& rand = utils::Randomize::instance();
            partitioned_hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
                partitioned_hypergraph.setOnlyNodePart(hn, rand.getRandomInt(0, k - 1, sched_getcpu()));
            });
        }

    } // namespace ds
} // namespace mt_kahypar


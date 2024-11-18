#pragma once

namespace mt_kahypar::dyn {

    // incoming change (nodes, edges, pins)
    struct Change {
        std::vector <HypernodeID> added_nodes{};
        std::vector <HyperedgeID> added_edges{};
        std::vector <HypernodeID> added_pins{};

        std::vector <HypernodeID> removed_nodes{};
        std::vector <HyperedgeID> removed_edges{};
        std::vector <HypernodeID> removed_pins{};
    };
}
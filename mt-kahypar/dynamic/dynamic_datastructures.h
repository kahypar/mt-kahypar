#pragma once

namespace mt_kahypar::dyn {

    // changing pin
    struct PinChange {
        HypernodeID node;
        HyperedgeID edge;
    };

    // incoming change (nodes, edges, pins)
    struct Change {
        std::vector <HypernodeID> added_nodes{};
        std::vector <HyperedgeID> added_edges{};
        std::vector <PinChange> added_pins{};

        std::vector <HypernodeID> removed_nodes{};
        std::vector <HyperedgeID> removed_edges{};
        std::vector <PinChange> removed_pins{};
    };
}
#pragma once

namespace mt_kahypar::dyn {

    // changing pin
    struct PinChange {
        HypernodeID node;
        HyperedgeID edge;

        // Overload << operator for PinChange
        friend std::ostream& operator<<(std::ostream& os, const PinChange& pin) {
          os << "(" << pin.node << "," << pin.edge << ")";
          return os;
        }
    };

    // incoming change (nodes, edges, pins)
    struct Change {
        std::vector <HypernodeID> added_nodes{};
        std::vector <HyperedgeID> added_edges{};
        std::vector <PinChange> added_pins{};

        std::vector <HypernodeID> removed_nodes{};
        std::vector <HyperedgeID> removed_edges{};
        std::vector <PinChange> removed_pins{};

        void append(const Change& other) {
            added_nodes.insert(added_nodes.end(), other.added_nodes.begin(), other.added_nodes.end());
            added_edges.insert(added_edges.end(), other.added_edges.begin(), other.added_edges.end());
            added_pins.insert(added_pins.end(), other.added_pins.begin(), other.added_pins.end());
            removed_nodes.insert(removed_nodes.end(), other.removed_nodes.begin(), other.removed_nodes.end());
            removed_edges.insert(removed_edges.end(), other.removed_edges.begin(), other.removed_edges.end());
            removed_pins.insert(removed_pins.end(), other.removed_pins.begin(), other.removed_pins.end());
        }

        // Overload << operator for Change
        friend std::ostream& operator<<(std::ostream& os, const Change& change) {
          os << "Added Nodes: ";
          for (const auto& node : change.added_nodes) {
            os << node << " ";
          }
          os << "\nAdded Edges: ";
          for (const auto& edge : change.added_edges) {
            os << edge << " ";
          }
          os << "\nAdded Pins: ";
          for (const auto& pin : change.added_pins) {
            os << pin << " ";
          }

          os << "\nRemoved Nodes: ";
          for (const auto& node : change.removed_nodes) {
            os << node << " ";
          }
          os << "\nRemoved Edges: ";
          for (const auto& edge : change.removed_edges) {
            os << edge << " ";
          }
          os << "\nRemoved Pins: ";
          for (const auto& pin : change.removed_pins) {
            os << pin << " ";
          }

          return os;
        }
    };
}
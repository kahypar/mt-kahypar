/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "datastructure/flow_hypergraph.h"

#include "mt-kahypar/macros.h"

namespace mt_kahypar {

	class FlowHypergraphBuilder : public whfc::FlowHypergraph {
	public:
		using Base = whfc::FlowHypergraph;

		FlowHypergraphBuilder() :
      Base(),
      _finalized(false),
      _numPinsAtHyperedgeStart(0) {
			clear();
		}

		explicit FlowHypergraphBuilder(size_t num_nodes) :
      Base(),
      _finalized(false),
      _numPinsAtHyperedgeStart(0) {
			reinitialize(num_nodes);
		}

    // ####################### Sequential Construction #######################

		void addNode(const whfc::NodeWeight w) {
			nodes.back().weight = w;
			nodes.push_back({whfc::InHeIndex(0), whfc::NodeWeight(0)});
		}

		void startHyperedge(const whfc::Flow capacity) {
			finishHyperedge();	//finish last hyperedge
			hyperedges.back().capacity = capacity;	//exploit sentinel
			_numPinsAtHyperedgeStart = numPins();
			maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, capacity);
		}

		void addPin(const whfc::Node u) {
			assert(u < numNodes());
			pins.push_back({u, whfc::InHeIndex::Invalid()});
			nodes[u+1].first_out++;
		}

		size_t currentHyperedgeSize() const {
			return numPins() - _numPinsAtHyperedgeStart;
		}

		void removeCurrentHyperedge() {
			while (numPins() > _numPinsAtHyperedgeStart) {
				removeLastPin();
      }
		}

		void finalize();

    // ####################### Parallel Construction #######################

    void allocateNodes(const size_t num_nodes) {
      nodes.assign(num_nodes + 1, NodeData { whfc::InHeIndex(0), whfc::NodeWeight(0) });
    }

    void resizeNodes(const size_t num_nodes) {
      ASSERT(num_nodes <= numNodes());
      nodes.resize(num_nodes + 1);
    }

    void allocateHyperedgesAndPins(const size_t num_hyperedges,
                                   const size_t num_pins);

    void resizeHyperedgesAndPins(const size_t num_hyperedges,
                                 const size_t num_pins);

    void addPin(const whfc::Node u, const size_t pin_idx) {
      ASSERT(pin_idx < pins.size());
      ASSERT(static_cast<size_t>(u) < numNodes());
      pins[pin_idx].pin = u;
      __atomic_fetch_add(reinterpret_cast<whfc::InHeIndex::ValueType*>(
        &nodes[u + 1].first_out), 1, __ATOMIC_RELAXED);
    }

    void finishHyperedge(const whfc::Hyperedge he, const whfc::Flow capacity,
                         const size_t pin_start_idx, const size_t pin_end_idx);

    void finalizeParallel();

    // ####################### Common Functions #######################

		void clear();

		void reinitialize(size_t num_nodes) {
			clear();
			nodes.resize(num_nodes + 1);
		}

		void shrink_to_fit() {
			nodes.shrink_to_fit();
			hyperedges.shrink_to_fit();
			pins.shrink_to_fit();
			incident_hyperedges.shrink_to_fit();
			pins_sending_flow.shrink_to_fit();
			pins_receiving_flow.shrink_to_fit();
		}

	private:

    // ####################### Sequential Construction #######################

		void removeLastPin() {
			nodes[ pins.back().pin + 1 ].first_out--;
			pins.pop_back();
		}

		bool finishHyperedge();

		bool _finalized;
		size_t _numPinsAtHyperedgeStart;

    // ####################### Parallel Construction #######################

    bool verifyParallelConstructedHypergraph();
	};
}
/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {
namespace io {

static inline void printBanner(const Context& context) {
  if (!context.partition.quiet_mode) {

    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
    LOG << R"(+         __  __ _______       _  __     _    _       _____                   +)";
    LOG << R"(+        |  \/  |__   __|     | |/ /    | |  | |     |  __ \                  +)";
    LOG << R"(+        | \  / |  | |  ____  | ' / __ _| |__| |_   _| |__) |_ _ _ __         +)";
    LOG << R"(+        | |\/| |  | | |____| |  < / _` |  __  | | | |  ___/ _` | '__|        +)";
    LOG << R"(+        | |  | |  | |        | . \ (_| | |  | | |_| | |  | (_| | |           +)";
    LOG << R"(+        |_|  |_|  |_|        |_|\_\__,_|_|  |_|\__, |_|   \__,_|_|           +)";
    LOG << R"(+                                                __/ |                        +)";
    LOG << R"(+                                               |___/                         +)";
    LOG << R"(+          Karlsruhe Shared Memory Hypergraph Partitioning Framework          +)";
    LOG << R"(+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++)";
  }
}

static inline void printInputInformation(const Context& context/*, const Hypergraph& hypergraph*/) {
  if (context.type == ContextType::main && !context.partition.quiet_mode) {
    LOG << context;
    /*if (context.partition.verbose_output) {
      LOG << "\n********************************************************************************";
      LOG << "*                                    Input                                     *";
      LOG << "********************************************************************************";
      io::printHypergraphInfo(hypergraph, context.partition.graph_filename.substr(
                                context.partition.graph_filename.find_last_of('/') + 1));
    }*/
  }
}

} // namespace io
} // namespace mt_kahypar

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

#include "kahypar/macros.h"

#define HEAVY_ASSERT_(TYPE, N) HEAVY_ ## TYPE ## _ASSERT_ ## N
#define HEAVY_ASSERT_EVAL(TYPE, N) HEAVY_ASSERT_(TYPE, N)

#ifdef KAHYPAR_USE_HEAVY_PREPROCESSING_ASSERTIONS
  #define HEAVY_PREPROCESSING_ASSERT_1(cond) ASSERT(cond)
  #define HEAVY_PREPROCESSING_ASSERT_2(cond, msg) ASSERT(cond, msg)
#else
  #define HEAVY_PREPROCESSING_ASSERT_1(cond)
  #define HEAVY_PREPROCESSING_ASSERT_2(cond, msg)
#endif

#ifdef KAHYPAR_USE_HEAVY_COARSENING_ASSERTIONS
  #define HEAVY_COARSENING_ASSERT_1(cond) ASSERT(cond)
  #define HEAVY_COARSENING_ASSERT_2(cond, msg) ASSERT(cond, msg)
#else
  #define HEAVY_COARSENING_ASSERT_1(cond)
  #define HEAVY_COARSENING_ASSERT_2(cond, msg)
#endif

#ifdef KAHYPAR_USE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS
  #define HEAVY_INITIAL_PARTITIONING_ASSERT_1(cond) ASSERT(cond)
  #define HEAVY_INITIAL_PARTITIONING_ASSERT_2(cond, msg) ASSERT(cond, msg)
#else
  #define HEAVY_INITIAL_PARTITIONING_ASSERT_1(cond)
  #define HEAVY_INITIAL_PARTITIONING_ASSERT_2(cond, msg)
#endif

#ifdef KAHYPAR_USE_HEAVY_REFINEMENT_ASSERTIONS
  #define HEAVY_REFINEMENT_ASSERT_1(cond) ASSERT(cond)
  #define HEAVY_REFINEMENT_ASSERT_2(cond, msg) ASSERT(cond, msg)
#else
  #define HEAVY_REFINEMENT_ASSERT_1(cond)
  #define HEAVY_REFINEMENT_ASSERT_2(cond, msg)
#endif

// Heavy assertions are assertions which increase the complexity of the scope
// which they are executed in by an polynomial factor. In debug mode you are often only
// interested in certain phase of the multilevel paradigm. However, when enabling all assertions
// it can take a while to reach the point which you are really interested in, because heavy assertions
// radicaly downgrade the performance of the application. Therefore such assertions such be packed
// in a heavy assertion macro. One can enable all heavy assertions of a certain phase via cmake flag
// KAHYPAR_HEAVY_ASSERTION_TYPE.
#define HEAVY_PREPROCESSING_ASSERT(...) EXPAND(HEAVY_ASSERT_EVAL(PREPROCESSING, EXPAND(NARG(__VA_ARGS__)))(__VA_ARGS__))
#define HEAVY_COARSENING_ASSERT(...) EXPAND(HEAVY_ASSERT_EVAL(COARSENING, EXPAND(NARG(__VA_ARGS__)))(__VA_ARGS__))
#define HEAVY_INITIAL_PARTITIONING_ASSERT(...) EXPAND(HEAVY_ASSERT_EVAL(INITIAL_PARTITIONING, EXPAND(NARG(__VA_ARGS__)))(__VA_ARGS__))
#define HEAVY_REFINEMENT_ASSERT(...) EXPAND(HEAVY_ASSERT_EVAL(REFINEMENT, EXPAND(NARG(__VA_ARGS__)))(__VA_ARGS__))
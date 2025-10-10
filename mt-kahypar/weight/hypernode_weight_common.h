/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2025 Nikolai Maas <nikolai.maas@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include "mt-kahypar/weight/allocated_hypernode_weight.h"
#include "mt-kahypar/weight/hypernode_weight_array.h"
#include "mt-kahypar/weight/hypernode_weight_base.h"
#include "mt-kahypar/weight/hypernode_weight_operators.h"

namespace mt_kahypar {

using mt_kahypar::weight::Dimension;

using mt_kahypar::weight::HNWeightScalar;
using mt_kahypar::weight::HNWeightRef;
using mt_kahypar::weight::HNWeightConstRef;
using mt_kahypar::weight::HNWeightAtomicCRef;
using mt_kahypar::weight::AllocatedHNWeight;

using mt_kahypar::weight::HypernodeWeightArray;

}  // namespace mt_kahypar

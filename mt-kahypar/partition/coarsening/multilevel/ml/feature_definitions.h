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

#include <array>

#include "kahypar-resources/meta/typelist.h"

#include "mt-kahypar/partition/coarsening/multilevel/ml/feature_mappers.h"

namespace mt_kahypar {
namespace features {
  using OrderedFeatures = kahypar::meta::Typelist<
                            comm_1_equal,
                            N1<to_n1_edges>,
                            modularity_1,
                            modularity_2,
                            N1<degree_quantile>,
                            N1<modularity>,
                            N0<min_contracted_degree>,
                            N0<to_n2_edges>,
                            modularity_0,
                            strawman_similarity,
                            N0<degree>,
                            N0<modularity>,
                            N0<q3_degree>,
                            N1<med_degree>,
                            N1<min_contracted_degree>,
                            N1<q1_degree>,
                            N0<med_degree>,
                            N0<skew_degree>,
                            N0<to_n1_edges>,
                            N1<q3_degree>,
                            N1<degree>,
                            N0<skew_degree>,
                            N1<skew_degree>,
                            N1<to_n2_edges>,
                            N0<sd_degree>,
                            N0<entropy_degree>,
                            N1<sd_degree>,
                            N1<avg_degree>,
                            N0<avg_degree>,
                            N0<degree_quantile>,
                            N0<q1_degree>,
                            N1<entropy_degree>>;
}

constexpr std::array<float, 32> MEANS = {
  5.587499999999999689e-01,
  1.340532353285714271e+04,
  4.246738004000000788e-01,
  4.083488801571429905e-01,
  7.678186083390997618e-01,
  4.597649769698764627e-02,
  1.708031885283828188e+01,
  8.120931144371428527e+05,
  4.203864180000000950e-01,
  9.032276219906050507e-02,
  1.852811740000000100e+03,
  6.640964825313666609e-02,
  5.927067171428570873e+02,
  7.735028757142856648e+02,
  1.686635149178571069e+01,
  1.994819099999999992e+02,
  2.062955028571428500e+02,
  6.724084657185714775e+01,
  1.664587893800000020e+05,
  1.479016372857142869e+03,
  1.694778171428571341e+02,
  5.960622753630457460e+00,
  2.201848402172734254e+00,
  2.706347131714285933e+05,
  1.149621365775856930e+03,
  9.024074388528566804e-01,
  1.988655032879999681e+03,
  1.454631654707142616e+03,
  5.631103718657142281e+02,
  9.444221769571425895e-01,
  7.066764142857142872e+01,
  8.986638262857142845e-01
};

constexpr std::array<float, 32> STDEVS = {
  4.965364412608606060e-01,
  8.588495500416113646e+04,
  2.672166113636331719e-01,
  2.549913137143254294e-01,
  2.426047473662316656e-01,
  1.178635320943709347e-01,
  1.551317710620404178e+01,
  1.967225653240529587e+06,
  2.625392904069992484e-01,
  2.798784014819127863e-01,
  6.938493833843884204e+03,
  1.464652059888615632e-01,
  1.175256826385265413e+03,
  3.615872713853701043e+03,
  1.654281096095436965e+01,
  1.111696943389394846e+03,
  3.757207162862234782e+02,
  1.045473265702798358e+02,
  7.292372495000237832e+05,
  4.714866372956524174e+03,
  4.719927170436072288e+02,
  1.072348905896073212e+01,
  1.955097228806399245e+00,
  8.446655574265166651e+05,
  2.087784866213510668e+03,
  1.647383706372066159e-01,
  4.529102928722440083e+03,
  4.238461564352129244e+03,
  9.758818421947960360e+02,
  1.102981501712899520e-01,
  1.228830924871626706e+02,
  2.766351295364234963e-01
};

}  // namespace mt_kahypar

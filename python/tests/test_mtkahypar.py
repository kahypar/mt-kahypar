#******************************************************************************
# MIT License
#
# This file is part of Mt-KaHyPar.
#
# Copyright (C) 2022 Tobias Heuer <tobias.heuer@kit.edu>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#*****************************************************************************/

import unittest
import os

import mtkahypar as mtkahypar

mydir = os.path.dirname(os.path.realpath(__file__))

class MainTest(unittest.TestCase):

  def test_hello_world(self):
    context = mtkahypar.Context()
    context.loadConfigurationFile("../../config/default_flow_preset.ini")
    context.setPartitioningParameters(8, 0.01, mtkahypar.Objective.CUT, 41)
    context.setIndividualBlockWeights([20, 30, 40, 20, 10, 30, 40, 50])
    context.setNumberOfVCycles(3)
    context.enableLogging(True)
    context.outputConfiguration()

if __name__ == '__main__':
  unittest.main()

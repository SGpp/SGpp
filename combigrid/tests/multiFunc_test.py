# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import sys, os
import pysgpp

def f(x):
    return x[0]*x[1]

##
# Class with unittests for multiFunc from combigrid module 
class TestMultifunc(unittest.TestCase):

    ## Set up the variables
    def setUp(self,):
        self.function = pysgpp.multiFunc(f)
        self.evalValue = pysgpp.DataVector(2)
        self.evalValue.set(0,2.2)
        self.evalValue.set(1,3)
        self.result = f(self.evalValue)

    ##
    # Tests the functions @link python.learner.formatter.GridFormatter.GridFormatter.serialize() GridFormatter.serialize() @endlink
    # and  @link python.learner.formatter.GridFormatter.GridFormatter.deserializeFromFile() GridFormatter.deserializeFromFile() @endlink
    def testEval(self,):
        res = self.function(self.evalValue)
        res_call = self.function.call(self.evalValue)
        self.assertTrue(res == self.result, "MultiFunction operator() failed")
        self.assertTrue(res_call == self.result, "MultiFunction operator call failed")


if __name__=="__main__":
    unittest.main()
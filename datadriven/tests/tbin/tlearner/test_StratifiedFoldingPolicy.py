###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Plueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Valeriy Khakhutskyy (khakhutv@in.tum.de)####################################################################

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.learner.folding.StratifiedFoldingPolicy import StratifiedFoldingPolicy
from bin.data.DataContainer import DataContainer
from pysgpp import DataVector, DataMatrix


##
# @package tests.tbin.test_StratifiedFoldingPolicy
# Contains class test_StratifiedFoldingPolicy::TestStratifiedFoldingPolicy with unittests for @link bin.learner.folding.StratifiedFoldingPolicy.StratifiedFoldingPolicy StratifiedFoldingPolicy @endlink

##
# Class with unittests for @link bin.learner.folding.StratifiedFoldingPolicy.StratifiedFoldingPolicy StratifiedFoldingPolicy @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.folding.StratifiedFoldingPolicy.StratifiedFoldingPolicy StratifiedFoldingPolicy @endlink
class TestStratifiedFoldingPolicy(unittest.TestCase):

    
    ## Set up the variables
    def setUp(self):
        self.size = 9
        self.level = 4
        points = DataMatrix(self.size, 1)
        values = DataVector(self.size)
        for i in xrange(self.size):
            points.set(i, 0, i)
            values[i] = -1 if i < self.size/2 else 1
        self.dataContainer = DataContainer(points=points, values=values)
        self.policy = StratifiedFoldingPolicy(self.dataContainer, self.level)
    
    
    ##
    # Tests the function @link bin.learner.folding.FoldingPolicy.FoldingPolicy.next() StratifiedFoldingPolicy.next() @endlink    
    def testNext(self):
        validationCorrectData = [[4,0],[5,1], [6,2], [7,8,3]]
        self.assertEqual(self.level, len(self.policy.dataFold))
        i = 0
        for l in self.policy:
            trainPoints = l.getTrainDataset().getPoints()
            validationPoints = l.getTestDataset().getPoints()
            validCorrectD = set(validationCorrectData[i])
            trainCorrectD = set(range(self.size)) - validCorrectD
            
            self.assertEqual(trainPoints.getNrows(), len(trainCorrectD))
            self.assertEqual(validationPoints.getSize(), len(validCorrectD))
            for k in xrange(trainPoints.getNrows()):
                self.assertTrue(trainPoints.get(k,0)  in trainCorrectD)
            for k in xrange(validationPoints.getSize()):
                self.assertTrue(validationPoints.get(k,0) in validCorrectD)
                
            i += 1
           
        
        
        
if __name__=="__main__":
    unittest.main() 
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.learner.folding import FilesFoldingPolicy
from pysgpp.extensions.datadriven.data.ARFFAdapter import ARFFAdapter


##
# @package tests.tbin.test_Classifier
# Contains class test_Classifier::TestClassifier with unittests for @link
# python.learner.folding.FilesFoldingPolicy.FilesFoldingPolicy FilesFoldingPolicy @endlink

##
# Class with unittests for @link
# python.learner.folding.FilesFoldingPolicy.FilesFoldingPolicy FilesFoldingPolicy @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.learner.folding.FilesFoldingPolicy.FilesFoldingPolicy FilesFoldingPolicy @endlink
class TestFilesFoldingPolicy(unittest.TestCase):

    ## Set up the variables
    def setUp(self):
        files =  ['/datasets/foldf_fold0.arff', '/datasets/foldf_fold1.arff', '/datasets/foldf_fold2.arff']
        datasets = []
        fileCounter = 0
        self.dataContainer = ARFFAdapter(pathlocal + files[fileCounter]).loadData("train" + str(fileCounter))
        for fname in files[1:]:
            fileCounter += 1
            self.dataContainer = self.dataContainer.combine(ARFFAdapter(pathlocal + fname).loadData("train" + str(fileCounter)))

        self.policy = FilesFoldingPolicy(self.dataContainer)

        self.points = list(range(9))
        self.values = list(range(9))
        self.values.reverse()


    ##
    # Tests the function @link python.learner.folding.FoldingPolicy.FoldingPolicy.next() FilesFoldingPolicy.next() @endlink
    def testNext(self):
#        validationCorrectData = [[4,0],[5,1], [6,2], [7,8,3]]
#        self.assertEqual(self.level, len(self.policy.dataFold))
        step = 0
        for l in self.policy:
            points = l.getPoints()
            testPoints = self.points[:step*3] + self.points[step*3+3:]
            values = l.getValues()
            testValues = self.values[:step*3] + self.values[step*3+3:]
            self.assertEqual(points.getNrows(), len(testPoints))
            self.assertEqual(len(values), len(testValues))

            for i in range(points.getSize()):
                self.assertEqual(points.get(i,0), testPoints[i])
                self.assertEqual(values[i], testValues[i])
            step += 1




if __name__=="__main__":
    unittest.main()

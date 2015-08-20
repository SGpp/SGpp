# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.learner.LearnerBuilder import LearnerBuilder
from pysgpp.extensions.datadriven.controller import InfoToScreen
from pysgpp.extensions.datadriven.controller import InfoToFile


##
# @package tests.tbin.test_LearnerBuilder
# Contains class test_LearnerBuilder::TestLearnerBuilder with unittests for @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder LearnerBuilder @endlink

##
# Class with unittests for @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder LearnerBuilder @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder LearnerBuilder @endlink
class TestLearnerBuilder(unittest.TestCase):


    ## Set up the variables and
    # test the functions @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder.buildClassifier() LearnerBuilder.buildClassifier() @endlink
    # @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder.withTrainingDataFromARFFFile() LearnerBuilder.withTrainingDataFromARFFFile() @endlink
    # @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder.withStopPolicy() LearnerBuilder.withStopPolicy() @endlink
    # @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder.withCGSolver() LearnerBuilder.withCGSolver() @endlink
    # @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder.withProgressPresenter() LearnerBuilder.withProgressPresenter() @endlink
    # @link python.pysgpp_datadriven.learner.LearnerBuilder.LearnerBuilder.andGetResult() LearnerBuilder.andGetResult() @endlink
    def setUp(self):
        self.builder = LearnerBuilder()
        self.classifier = self.builder.buildClassifier().withTrainingDataFromARFFFile(pathlocal + "/datasets/classifier.train.arff")\
                    .withTestingDataFromARFFFile(pathlocal + "/datasets/classifier.test.arff").withGrid().withLevel(2).withSpecification().withLambda(0.00001).withAdaptPoints(2)\
                    .withStopPolicy()\
                    .withAdaptiveItarationLimit(1).withCGSolver().withProgressPresenter(InfoToFile(pathlocal + "/presentor.test"))\
                    .andGetResult()
#        level = 2
#        dim = 2
#        l = 0.00001
#        self.classifier = Classifier()
#        dataContainer = ARFFAdapter(pathlocal + "/datasets/classifier.train.arff").loadData()
#        self.classifier.setDataContainer(dataContainer)
#        foldingPolicy = FoldingPolicy(dataContainer)
#        self.classifier.setFoldingPolicy(foldingPolicy)
#        grid = Grid.createLinearGrid(dim)
#        storage = grid.createGridGenerator()
#        storage.regular(level)
#        self.classifier.setGrid(grid)
#        self.classifier.setLearnedKnowledge(LearnedKnowledge(None))
#        spec = TrainingSpecification()
#        spec.setL(l)
#        self.classifier.setSpecification(spec)
#        stopPolicy = TrainingStopPolicy()
#        stopPolicy.setAdaptiveIterationLimit(0)
#        self.classifier.setStopPolicy(stopPolicy)
#        self.classifier.setSolver(CGSolver())


#    def testScreenEventPresentor(self,):
#        self.classifier.learnDataWithTest()



if __name__=="__main__":
    unittest.main()

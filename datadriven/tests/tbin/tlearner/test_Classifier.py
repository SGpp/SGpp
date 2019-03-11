# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.data.ARFFAdapter import ARFFAdapter
from pysgpp import *
from pysgpp.extensions.datadriven.learner.solver import CGSolver
from pysgpp.extensions.datadriven.learner import Learner, LearnedKnowledge, TrainingSpecification, TrainingStopPolicy, Classifier
from pysgpp.extensions.datadriven.learner.folding import SequentialFoldingPolicy, FoldingPolicy
from pysgpp.extensions.datadriven.controller.InfoToScreen import InfoToScreen


##
# @package tests.tbin.test_Classifier
# Contains class test_Classifier::TestClassifier with unittests for @link python.learner.Classifier.Classifier Classifier @endlink

##
# Class with unittests for @link python.learner.Classifier.Classifier Classifier @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.learner.Classifier.Classifier Classifier @endlink
class TestClassifier(unittest.TestCase):

    ## Set up the variables
    def setUp(self):
        level = 2
        dim = 2
        l = 0.00001
        self.classifier = Classifier()
        dataContainer = ARFFAdapter(pathlocal + "/datasets/classifier.train.arff").loadData()
        self.classifier.setDataContainer(dataContainer)
        foldingPolicy = FoldingPolicy(dataContainer)
        self.classifier.setFoldingPolicy(foldingPolicy)
        grid = Grid.createLinearGrid(dim)
        storage = grid.getGenerator()
        storage.regular(level)
        self.classifier.setGrid(grid)
        self.classifier.setLearnedKnowledge(LearnedKnowledge())
        spec = TrainingSpecification()
        spec.setL(l)
        spec.setCOperator(createOperationLaplace(grid))
        spec.setCOperatorType('laplace')
        self.classifier.setSpecification(spec)
        stopPolicy = TrainingStopPolicy()
        stopPolicy.setAdaptiveIterationLimit(0)
        self.classifier.setStopPolicy(stopPolicy)
        solver = CGSolver()
        #solver.attachEventController(InfoToScreen())
        solver.setImax(500)
        solver.setReuse(True)
        self.classifier.setSolver(solver)

#     ##
#     # Tests the function @link python.learner.Learner.Learner.learnDataWithTest() Classifier.learnDataWithTest() @endlink
#     def testLearnDataWithTest(self,):
#         correct = [-0.33360635579319858346,
#                    0.67890792146517364714,
#                    17.37781054927400248289,
#                    19.86707480839170614217,
#                    -1.61960456623131343612,
#                    -7.70442957659598182119,
#                    -22.14900166819601423640,
#                    8.92745373469135827804,
#                    -11.12477960213342775830
#                    ]
#         testDataset = ARFFAdapter(pathlocal + "/datasets/classifier.test.arff").loadData("test")
#         self.classifier.setDataContainer(self.classifier.dataContainer.combine(testDataset))
#         self.classifier.stopPolicy.setAdaptiveIterationLimit(1)
#         self.classifier.specification.setAdaptPoints(1)
#         alpha = self.classifier.learnDataWithTest(self.classifier.dataContainer.combine(testDataset))
#         self.assertEqual(len(alpha), len(correct))
#         for i in xrange(len(correct)):
#             self.assertAlmostEqual(alpha[i], correct[i], 3)

    ##
    # Tests the function @link python.learner.Learner.Learner.applyData() Classifier.applyData() @endlink
    def testApplyData(self,):
        correct = [0.253400605292, -0.25507958758, 0.0530555506998]
        points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
        self.classifier.learnData()
        data = DataMatrix(3,2)
        for i in range(3):
            temp = DataVector(2)
            temp[0] = points[i][0]
            temp[1] = points[i][1]
            data.setRow(i, temp)

        val = self.classifier.applyData(data)
        places = 7

        self.assertEqual(len(val), len(correct))
        for i in range(len(correct)):
            self.assertAlmostEqual(val[i], correct[i], places=places)

    ##
    # Tests the function @link python.learner.Learner.Learner.learnData() Classifier.learnData() @endlink
    def testLearnData(self):
        correct = [-0.03105750236900508068,
                   -0.61865507797281660274,
                   0.64903026441541222802,
                   0.64903026441541267211,
                   -0.61865507797281593660]

        alpha = self.classifier.learnData()
        places = 7
        for i in range(len(alpha)):
            self.assertAlmostEqual(correct[i], alpha[i], places=places)

    ##
    # Tests the function @link python.learner.Learner.Learner.learnDataWithFolding() Classifier.learnDataWithFolding() @endlink
    def testLearnDataWithFolding(self,):
        correct = [0.6612903226, 0.1428571429,
                   0.5741935484, 0.9142857143,
                   0.6193548387, 0.5142857143,
                   0.5870967742, 0.7714285714,
                   0.6032258065, 0.5714285714,
                   0.6387096774, 0.4000000000,
                   0.5935483871, 0.7428571429,
                   0.6193548387, 0.5142857143,
                   0.5903225806, 0.7714285714,
                   0.6063492063, 0.5666666667]
        level = 2
        dim = 6
        l = 0.00001
        self.classifier = Classifier()
        dataContainer = ARFFAdapter(pathlocal + "/datasets/liver-disorders_normalized.arff").loadData()
        self.classifier.setDataContainer(dataContainer)
        foldingPolicy = SequentialFoldingPolicy(dataContainer, 10)
        self.classifier.setFoldingPolicy(foldingPolicy)
        grid = Grid.createLinearGrid(dim)
        storage = grid.getGenerator()
        storage.regular(level)
        self.classifier.setGrid(grid)
        self.classifier.setLearnedKnowledge(LearnedKnowledge())
        spec = TrainingSpecification()
        spec.setL(l)
        spec.setCOperator(createOperationLaplace(grid))
        spec.setCOperatorType('laplace')
        self.classifier.setSpecification(spec)
        stopPolicy = TrainingStopPolicy()
        stopPolicy.setAdaptiveIterationLimit(0)
        self.classifier.setStopPolicy(stopPolicy)
        self.classifier.setSolver(CGSolver())

        self.classifier.learnDataWithFolding()
        for i in range(10):
            self.assertAlmostEqual(correct[2*i], self.classifier.trainAccuracy[i])
            self.assertAlmostEqual(correct[2*i+1], self.classifier.testAccuracy[i])


if __name__=="__main__":
    unittest.main(verbosity=9)

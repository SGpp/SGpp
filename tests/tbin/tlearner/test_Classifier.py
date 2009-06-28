##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)               #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################


import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.data.ARFFAdapter import ARFFAdapter
from bin.pysgpp import *
from bin.learner import CGSolver, Learner, LearnedKnowledge, TrainingSpecification, TrainingStopPolicy, FoldingPolicy, Classifier
from bin.learner import SequentialFoldingPolicy

class TestClassifier(unittest.TestCase):
    
    classifier = None
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
        storage = grid.createGridGenerator()
        storage.regular(level)
        self.classifier.setGrid(grid)
        self.classifier.setLearnedKnowledge(LearnedKnowledge(None))
        spec = TrainingSpecification()
        spec.setL(l)
        self.classifier.setSpecification(spec)
        stopPolicy = TrainingStopPolicy()
        stopPolicy.setAdaptiveIterationLimit(0)
        self.classifier.setStopPolicy(stopPolicy)
        self.classifier.setSolver(CGSolver())
       # self.classifier.
    
    def testLearnDataWithTest(self,):
        correct = [-0.33360635579319858346,
                   0.67890792146517364714,
                   17.37781054927400248289,
                   19.86707480839170614217,
                   -1.61960456623131343612,
                   -7.70442957659598182119,
                   -22.14900166819601423640,
                   8.92745373469135827804,
                   -11.12477960213342775830
                   ]
        testDataset = ARFFAdapter(pathlocal + "/datasets/classifier.test.arff").loadData("test")
        self.classifier.setDataContainer(self.classifier.dataContainer.combine(testDataset))
        self.classifier.stopPolicy.setAdaptiveIterationLimit(1)
        self.classifier.specification.setAdaptPoints(1)
        alpha = self.classifier.learnDataWithTest(self.classifier.dataContainer.combine(testDataset))
        self.assertEqual(alpha.getSize(), len(correct))
        for i in xrange(len(correct)):
            self.assertAlmostEqual(alpha[i], correct[i], 3)

    def testApplyData(self,):
        correct = [0.253400605292, -0.25507958758, 0.0530555506998]
        points = [[0.5, 0.1], [0.3, 0.4], [0.9, 0.7]]
        self.classifier.learnData()
        data = DataVector(3,2)
        for i in xrange(3):
            temp = DataVector(2)
            temp[0] = points[i][0]
            temp[1] = points[i][1]
            data.setRow(i, temp)
        val = self.classifier.applyData(data)
        self.assertEqual(val.getSize(), len(correct))
        for i in xrange(len(correct)):
            self.assertAlmostEqual(val[i], correct[i])
    
#    def testCalcResult(self,):
#        self.fail("Not implemented")

    def testLearnData(self):
        correct = [-0.03105750236900508068, 
                   -0.61865507797281660274, 
                   0.64903026441541222802,
                   0.64903026441541267211,
                   -0.61865507797281593660]

        alpha = self.classifier.learnData()
        for i in xrange(alpha.getSize()):
            self.assertAlmostEqual(correct[i], alpha[i])

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
        storage = grid.createGridGenerator()
        storage.regular(level)
        self.classifier.setGrid(grid)
        self.classifier.setLearnedKnowledge(LearnedKnowledge(None))
        spec = TrainingSpecification()
        spec.setL(l)
        self.classifier.setSpecification(spec)
        stopPolicy = TrainingStopPolicy()
        stopPolicy.setAdaptiveIterationLimit(0)
        self.classifier.setStopPolicy(stopPolicy)
        self.classifier.setSolver(CGSolver())
        
        self.classifier.learnDataWithFolding()
        for i in xrange(10):
            self.assertAlmostEqual(correct[2*i], self.classifier.trainAccuracy[i])
            self.assertAlmostEqual(correct[2*i+1], self.classifier.testAccuracy[i])
        
        
#    def testEvalError(self, ):
#        self.fail("Not implemented")
    
#    def testAddToRefine(self,):
#        self.fail("Not implemented")
    
#    def testUpdateResults(self,):
#        self.fail("Not implemented")
        
if __name__=="__main__":
    unittest.main() 

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

from bin.learner.LearnerBuilder import LearnerBuilder
from bin.controller import InfoToScreen
from bin.controller import InfoToFile


class TestLearnerBuilder(unittest.TestCase):
    
    builder = None
    classifier = None
    
    def setUp(self):
        self.builder = LearnerBuilder()
        self.classifier = self.builder.buildClassifier().withTrainingDataFromARFFFile(pathlocal + "/datasets/classifier.train.arff")\
                    .withTestingDataFromARFFFile(pathlocal + "/datasets/classifier.test.arff").withGrid().withLevel(2).withSpecification().withLambda(0.00001).withAdaptPoints(2)\
                    .withStopPolicy()\
                    .withAdaptiveItarationLimit(1).withCGSolver().withProgressPresentor(InfoToFile(pathlocal + "/presentor.test"))\
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

    
    def testScreenEventPresentor(self,):
        self.classifier.learnDataWithTest()


        
if __name__=="__main__":
    unittest.main() 

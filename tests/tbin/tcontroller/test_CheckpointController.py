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
from tbin.tdata.test_ARFFAdapter import pathlocal

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.pysgpp import Grid, DataVector
from bin.learner import LearnedKnowledge
import bin.learner.LearnerBuilder as LearnerBuilder

import gzip

from bin.controller import CheckpointController

class TestCheckpointController(unittest.TestCase):

    def testSaveGrid(self):
        dim = 2
        level = 2
        grid = Grid.createLinearGrid(dim)
        generator = grid.createGridGenerator()
        generator.regular(level)
        
        controller = CheckpointController("savegrid")
        controller.setGrid(grid)
        controller.saveGrid(0)
        
        f = gzip.open(pathlocal + "/savegrid.0.grid.gz", "r")
        try:
            sampleString = f.read()
        finally:
            f.close()
            
        self.assertEqual(grid.serialize(), sampleString)

        
    def testLoadGrid(self):
        dim = 2
        level = 2
        grid = Grid.createLinearGrid(dim)
        generator = grid.createGridGenerator()
        generator.regular(level)
        
        controller = CheckpointController("sample", pathlocal)
        sampleGrid = controller.loadGrid(0)
        
        # check dimension and size
        self.assertEqual(dim, sampleGrid.getStorage().dim())
        self.assertEqual(grid.getStorage().size(), sampleGrid.getStorage().size())
        
        # if string representations are equal, then grids are equal
        self.assertEqual(grid.serialize(), sampleGrid.serialize())
        
    def testSaveLearnedKnowledge(self):
        testValues = [-0.0310651210442,
                      -0.618841896127,
                       0.649230972775,
                       0.649230972775,
                      -0.618841896127]
        alpha = DataVector(len(testValues))
        for i in xrange(len(testValues)):
            alpha[i] = testValues[i]
        
        learnedKnowledge = LearnedKnowledge()
        learnedKnowledge.update(alpha)
        
        controller = CheckpointController("saveknowledge", pathlocal)
        controller.setLearnedKnowledge(learnedKnowledge)
        controller.saveLearnedKnowledge(0)
        
        f = gzip.open(pathlocal + "/saveknowledge.0.arff.gz", "r")
        try:
            sampleLines = f.readlines()[3:8]
        finally:
            f.close()
            
        self.assertEqual(testValues, [float(i) for i in sampleLines])
        
        
    def testLoadLearnedKnowledge(self):
        controller = CheckpointController("sample", pathlocal)
        learnedKnowledge = controller.loadLearnedKnowledge(0)
        
        testValues = [-0.0310651210442,
                      -0.618841896127,
                       0.649230972775,
                       0.649230972775,
                      -0.618841896127]
        
        self.assertEqual(1, learnedKnowledge.getAlphas().getDim())
        self.assertEqual(len(testValues), learnedKnowledge.getAlphas().getSize())
                         
        for i in xrange(len(testValues)):
            self.assertAlmostEqual(testValues[i], learnedKnowledge.getAlphas()[i])
            
    def testSaveAllLoadAll(self):
        # test of two method is put together since it should test the capability 
        # to store and restore data accurately
        builder = LearnerBuilder()

        # as storing of grid and knowledge is covered with other tests, only the test of learner is relevant and combination is
        classifier = builder.buildClassifier()\
                     .withTrainingDataFromARFFFile(pathlocal + "/traindata.arff")\
                     .withGrid().withLevel(2)\
                     .withSpecification().withLambda(0.00001).withAdaptPoints(2)\
                     .withStopPolicy().withAdaptiveItarationLimit(1)\
                     .withCGSolver().withImax(500)\
                     .andGetResult()
        classifier.learnData()
        
        string1 = classifier.toString()
        #print "string1:  ",string1
        
        controller = CheckpointController("saveload", pathlocal)

        controller.setLearner(classifier)
        controller.saveAll(0)
        
        del controller
        
        controller = CheckpointController("saveload", pathlocal)
        newClassifier = controller.loadAll(0)
        
        # quick and dirty way to compare to objects - with their string representation, would work only
        # if toString() method works properly

        string2 = newClassifier.toString()
        
        import difflib
        for line in difflib.context_diff(string1.split("\n"), string2.split("\n")):
            print line
        #print "string2:  ",string2
        self.assertEqual(classifier.toString(), newClassifier.toString())


            
if __name__=="__main__":
    unittest.main()
        
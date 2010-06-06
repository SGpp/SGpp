##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2008 Dirk Plueger (pflueged@in.tum.de)                      #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
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

from bin.learner.formatter.LearnedKnowledgeFormatter import LearnedKnowledgeFormatter
from bin.learner.LearnedKnowledge import LearnedKnowledge
from bin.pysgpp import DataVector


##
# @package tests.tbin.test_LearnedKnowledgeFormatter
# Contains class test_LearnedKnowledgeFormatte::TestLearnedKnowledgeFormatter with unittests for @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink

##
# Class with unittests for @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink
class TestLearnedKnowledgeFormatter(unittest.TestCase):
    
    ## Set up the variables
    def setUp(self,):
        self.filename_load = pathlocal + "/datasets/load.alpha.arff"
        self.filename_save = pathlocal + "/datasets/save.alpha.arff"
        self.formatter = LearnedKnowledgeFormatter()

    
    ##
    # Tests the function @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter.deserializeFromFile() LearnedKnowledgeFormatter.deserializeFromFile() @endlink
    def testLoad(self,):
        alphas = self.formatter.deserializeFromFile(self.filename_load)
        self.assertEqual(alphas.getSize(), 10)
        self.assertEqual(alphas.getDim(), 1)
        a = 0.1
        row = DataVector(1)
        for i in xrange(10):
            alphas.getRow(i, row)
            self.assertAlmostEqual(row[0], a)
            a = a + 0.1
    
    
    ##
    # Tests the functions @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter.serializeToFile() LearnedKnowledgeFormatter.serializeToFile() @endlink    
    def testSave(self,):
        alphas = DataVector(10)
        a = 0.1
        for i in xrange(10):
            alphas[i] = a
            a = a + 0.1
        knowledge = LearnedKnowledge()
        knowledge.update(alphas)
        self.formatter.serializeToFile(knowledge.createMemento(), self.filename_save)
        
        alphas = self.formatter.deserializeFromFile(self.filename_save)
        self.assertEqual(alphas.getSize(), 10)
        self.assertEqual(alphas.getDim(), 1)
        a = 0.1
        row = DataVector(1)
        for i in xrange(10):
            alphas.getRow(i, row)
            self.assertAlmostEqual(row[0], a)
            a = a + 0.1
        
        
        
if __name__=="__main__":
    unittest.main() 
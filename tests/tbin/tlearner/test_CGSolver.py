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
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.learner.CGSolver import CGSolver
from bin.pysgpp import *

class TestCGSolver(unittest.TestCase):
    
    linSys = None
    bOp = None
    cOp = None
    data = None
    grid = None
    y = None
    solver = None
    
    def setUp(self):
#        data = DataVector(3,2)
#        t = DataVector(2)
#        t[0] = 0.25
#        t[1] = 0.25
#        data.setRow(0,t)
#        t[0] = 0.5
#        t[1] = 0.5
#        data.setRow(1,t)
#        t[0] = 0.75
#        t[1] = 0.5
#        data.setRow(2,t)
#        self.data = data
#        self.grid = Grid.createLinearGrid(2)
#        storage = self.grid.getStorage()
#        generator = self.grid.createGridGenerator()
#        generator.regular(2)
#        self.bOp = self.grid.createOperationB()
#        self.cOp = self.grid.createOperationLaplace()
#        self.y = DataVector(3)   
#        self.y[0] = 1
#        self.y[1] = 1
#        self.y[2] = 1
#        self.linSys = LinearSystem(data, self.y, self.bOp, self.cOp, self.grid, 10**(-5))
        
        self.solver = CGSolver()
        self.solver.setImax(400)
        self.solver.setAccuracy(0.000001)
        print self.solver.getAccuracy()
    
    def testSetAccuracy(self):
        self.solver.setAccuracy(0.000001)
        print self.solver.getAccuracy()
        
#    def testSolve(self):
#        calculated = self.solver.solve(self.linSys)
#        correct = [ 0.999967, 0.749828, 0.499951, 0.749828, -5.99811e-05 ]
#        self.assertEqual(len(calculated), len(correct))
#        for i in xrange(len(correct)):
#            self.assertAlmostEqual(calculated[i], correct[i], 5)
    
if __name__=="__main__":
    unittest.main() 

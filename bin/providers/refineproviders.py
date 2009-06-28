# This file is part of sgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#

from bin.pysgpp import *

class RefineProvider:
    gridDecorator = None
    def __init__(self, gridDecorator):
        self.gridDecorator = gridDecorator
        print "constructor called"
        print "size: %d" %self.gridDecorator.getSize()
        
    def reset(self):
        """Called before each learning step"""
        self.alpha = DataVector(self.gridDecorator.getSize())
        self.alpha.setAll(0.0)
        self.i = 0
        
    def refine(self, adapt_points):
        """Called after learning all folds"""
        self.alpha.mult(1.0/self.i)
        self.gridDecorator.grid.createGridGenerator().refine(SurplusRefinementFunctor(self.alpha, adapt_points))
        print "GridPoints: ", self.gridDecorator.getSize()
        
    
    def add(self, alpha, training):
        """Called once per fold"""
        self.alpha.add(alpha)
        self.i += 1
        
class SurplusRefineProvider(RefineProvider):

    def __init__(self, gridDecorator):
        RefineProvider.__init__(self, gridDecorator)

##
# @todo rework
class ErrorRefineProvider(RefineProvider):
    
    def add(self, alpha, training):
        """Called once per fold"""
        
        temp = DataVector(training[1].getSize())
        B = self.gridDecorator.grid.createOperationB()
        
        B.multTranspose(alpha, training[0], temp)
        
        temp.sub(training[1])
        temp.sqr()        
        
        error_temp = DataVector(alpha.getSize())
        B.mult(temp, training[0], error_temp)
        
        RefineProvider.add(error_temp, training)

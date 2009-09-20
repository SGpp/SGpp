##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2008 Dirk Plueger (pflueged@in.tum.de)                      #
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

## @package CGSolver
# @ingroup bin.learner
# @brief Abstract class for solving of linear equations
# @version $CURR$

from LinearSolver import LinearSolver, LinearSolverEvents
from bin.pysgpp import *
from bin.learner.LinearSolver import LinearSolver

class CGSolver(ConjugateGradients, LinearSolver):

    accuracy = 0.0001   #the relationship of the norm of end residual to the normal of initial residual
    imax = 400          #maximal number of iterations used in CG
    delta_0 = 0         #norm of initial residual
    delta_new = 0       #norm of current residual
    alpha = None        #result vector
    
    def __init__(self,):
        ConjugateGradients.__init__(self, self.imax, self.accuracy)
        LinearSolver.__init__(self)


    ##Sets the accuracy parameter
    #
    # @param accuracy: float value of accuracy parameter
    # @return: CG Solver itself
    def setAccuracy(self, accuracy):
        ConjugateGradients.setEpsilon(self, accuracy)
        return self
    
    
    def getAccuracy(self):
        return self.myEpsilon
    
    ##Sets the maximal number of iterations
    #
    # @param imax: integer limit of number of iterations
    # @return: SG Solver itself 
    def setImax(self, imax):
        ConjugateGradients.setMaxIterations(self, imax)
        return self
    
    
    def getImax(self):
        return self.nMaxIterations
        
        
    def starting(self):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.STARTING)
    
    
    def calcStarting(self, ):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.CALC_STARTING)
    
    
    def iterationComplete(self, ):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.ITERATION_COMPLETE)
    
    
    def complete(self, ):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.COMPLETE)

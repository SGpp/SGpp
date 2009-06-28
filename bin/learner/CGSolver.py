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


class CGSolver(LinearSolver):

    accuracy = 0.0001   #the relationaship of the norm of end residual to the normal of initial residual
    imax = 400          #maximal number of iterations used in CG
    delta_0 = 0         #norm of initial residual
    delta_new = 0       #norm of current residual
    alpha = None        #result vector


    ##Sets the accuracy parameter
    #
    # @param accuracy: float value of accuracy parameter
    # @return: CG Solver itself
    def setAccuracy(self, accuracy):
        self.accuracy = accuracy
        return self
    
    
    ##Sets the maximal number of iterations
    #
    # @param imax: integer limit of number of iterations
    # @return: SG Solver itself 
    def setImax(self, imax):
        self.imax = imax
    
    
    ##Solve linear system with CG Method
    #
    # @param linearSystem: LinearSystem to solve
    # @return: DataVector with solution 
    def solve(self, linearSystem): #taken from new_cg() in painlesscg.py
        self.notifyEventControllers(LinearSolverEvents.STARTING)

        epsilon2 = self.accuracy*self.accuracy
        
        self.alpha = DataVector(linearSystem.getNumGridPoints())
        self.alpha.setAll(0.0)
        
        b = DataVector(linearSystem.getNumGridPoints())
        linearSystem.getRightHandSide(b)
        
        self.iteration = 1
        temp = DataVector(self.alpha.getSize())
        q = DataVector(self.alpha.getSize())
        
        self.delta_0 = 0.0
        
        # calculate residuum
#        if reuse:
#            q.setAll(0)
#            linearSystem.apply(q, temp)
#            r = DataVector(b)
#            r.sub(temp)
#            self.delta_0 = r.dotProduct(r)*epsilon2
#        else:
        self.alpha.setAll(0) #@fixme: alpha reusing so far not supported
    
        linearSystem.apply(self.alpha, temp)
        r = DataVector(b)
        r.sub(temp)
    
        # delta
        d = DataVector(r)
        
        delta_old = 0.0
        self.delta_new = r.dotProduct(r)
    
#        if not reuse:
        self.delta_0 = self.delta_new*epsilon2 #@fixme: alpha reuse so far not supported
        
        self.notifyEventControllers(LinearSolverEvents.CALC_STARTING)
    
        while (self.iteration < self.imax+1) and (self.delta_new > self.delta_0):
            # q = A*d
            linearSystem.apply(d, q)
            # a = d_new / d.q
            a = self.delta_new/d.dotProduct(q)
    
            # x = x + a*d
            self.alpha.axpy(a, d)
    
            if self.iteration % 50 == 0:
            # r = b - A*x
                linearSystem.apply(self.alpha, temp)
                r.copyFrom(b)
                r.sub(temp)
            else:
                # r = r - a*q
                r.axpy(-a, q)
            
            self.notifyEventControllers(LinearSolverEvents.ITERATION_COMPLETE)
            
            delta_old = self.delta_new
            self.delta_new = r.dotProduct(r)
            beta = self.delta_new/delta_old
            
            d.mult(beta)
            d.add(r)
            
            self.iteration += 1
            
        self.notifyEventControllers(LinearSolverEvents.COMPLETE)
        return self.alpha



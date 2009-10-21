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


from LinearSolver import LinearSolver, LinearSolverEvents
from bin.pysgpp import *
from bin.learner.LinearSolver import LinearSolver
import types

## This is a <a href="http://en.wikipedia.org/wiki/Decorator_pattern" target="new">decorator</a> for 
# @link sg::solver::sle::ConjugateGradients ConjugateGradients@endlink class.
# The ConjugateGradients solver is enhanced with methods of concrete subject of <a href="http://en.wikipedia.org/wiki/Observer_pattern" target="new">the observer design pattern</a>
# described in @link bin.learner.LinearSolver.LinearSolver LinearSolver@endlink and function for serialization
# end deserialization.
#@todo: (khakhutv) rename set/getAccuracy and set/getImax for consistency with ConjugateGradients (low)
#@todo (khakhutv) currently there is always two parameters for delta/residuum accuracy/myEpsilon imax/nMaxIterations
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
        
        
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object. 
    def toString(self):
        serializationString = "'module' : '" + self.__module__ + "',\n"
        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            
            #lists, integers, floats, dictionaries can serialized with str()
            if type(attrValue) in [types.ListType, types.IntType, 
                             types.FloatType, types.DictType] and attrName.find("__") != 0: 
                serializationString += "'" + attrName + "'" + " : " + str(attrValue) + ",\n"
            # serialize strings with quotes    
            elif type(attrValue) == types.StringType and attrName.find("__") != 0:
                serializationString += "'" + attrName + "'" + " ': " + attrValue + "',\n"
                
        serializationString = "{" + serializationString.rstrip(",\n") + "}\n"
        return serializationString    
    
    
    # Restores the CGSolver object from the json object with attributes.
    #
    # @param jsonObject A json object.
    # @return The restored SGSolver object.
    @classmethod
    def fromJson(cls, jsonObject):
        cg = CGSolver()
        if jsonObject.has_key('accuracy'):
            cg.setAccuracy( jsonObject['accuracy'] )
        if jsonObject.has_key('imax'):
            cg.setImax( jsonObject['imax'] )
        if jsonObject.has_key('delta_0'):
            cg.delta_0 = jsonObject['delta_0']
        if jsonObject.has_key('delta_new'):
            cg.delta_new = jsonObject['delta_new']
        if jsonObject.has_key('nIterations'):
            cg.nIterations = jsonObject['nIterations']
        if jsonObject.has_key('nMaxIterations'):
            cg.nMaxIterations = jsonObject['nMaxIterations']
        if jsonObject.has_key('residuum'):
            cg.residuum = float(jsonObject['residuum'])
        return cg
         
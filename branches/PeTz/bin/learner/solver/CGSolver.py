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
import types

## This is a <a href="http://en.wikipedia.org/wiki/Decorator_pattern" target="new">decorator</a> for sg::ConjugateGradients class.
# The ConjugateGradients solver is enhanced with methods of concrete subject of <a href="http://en.wikipedia.org/wiki/Observer_pattern" target="new">the observer design pattern</a>
# described in @link bin.learner.solver.LinearSolver LinearSolver@endlink and function for serialization
# end deserialization.
#
# In order to combine high performance of C++ code and flexibility of Subscription pattern
# the <a href="http://en.wikipedia.org/wiki/Template_pattern" target="new">Template design pattern</a>
# was used in this class. So the CG algorithm itself is implemented in C++ class
# ConjugateGradients, where template methods starting(), calcStarting(),
# iterationComplete(), and complete() are defined and called in different phases
# of the CG algorithm. This methods are overridden by CGSolver to rise the corresponding events by
# event subscribers.
class CGSolver(ConjugateGradients, LinearSolver):
    
    ##the relationship of the norm of end residual to the normal of initial residual
    DEFAULT_ACCURACY = 0.0001  
    
    ##maximal number of iterations used in CG
    DEFAULT_IMAX = 400          
    
    ##norm of initial residual
    #(khakhutv) I'm not sure if I should delete the start residuum, as it's not used anywhere
#    delta_0 = 0         
    
    ##norm of current residual
    #(khakhutv) I'm not sure if I should delete the start residuum, as there is ConjugateGradients.getResiduum() 
#    delta_new = 0    
    
    ##result vector   
    alpha = None  
    
        
    
    ##whether the old alpha vector should be reused
    __reuse = False  
    
    
    ## Constructor
    def __init__(self,):
        ConjugateGradients.__init__(self, self.DEFAULT_IMAX, self.DEFAULT_ACCURACY)
        LinearSolver.__init__(self)
        
        ##Maximal accuracy. If the norm of the residuum falls below max_threshold, stop the CG iterations. Default value: -1
        self.max_threshold = -1

    ## Returns True if the old alpha vector should be reused
    # @return True if the old alpha vector should be reused
    def getReuse(self):
        return self.__reuse

    ## Defines whether the old alpha vector should be reused
    # @param value True if the old alpha vector should be reused
    def setReuse(self, value):
        self.__reuse = value

    
    ## Returns the @link bin::learner::solver::CGSolver::CGSolver.max_threshold threshold@endlink parameter
    # @return the @link bin::learner::solver::CGSolver::CGSolver.max_threshold threshold@endlink parameter 
    def getThreshold(self):
        return self.max_threshold
    
    
    ##Sets the @link bin::learner::solver::CGSolver::CGSolver.max_threshold threshold@endlink parameter
    #
    # @param threshold: float value of @link bin::learner::solver::CGSolver::CGSolver.max_threshold threshold@endlink parameter
    # @return: CG Solver itself
    def setThreshold(self, threshold):
        self.max_threshold = threshold

    ##Sets the accuracy parameter
    #
    # @param accuracy: float value of DEFAULT_ACCURACY parameter
    # @return: CG Solver itself
    def setEpsilon(self, accuracy):
        ConjugateGradients.setEpsilon(self, accuracy)
        return self
    
    ## Return the accuracy for CG divergence criterion
    # @return the accuracy for CG divergence criterion
    def getEpsilon(self):
        return self.myEpsilon
    
    ##Sets the maximal number of iterations
    #
    # @param imax: integer limit of number of iterations
    # @return: SG Solver itself 
    def setImax(self, imax):
        #ConjugateGradients.setMaxIterations(self, imax)
        self.setMaxIterations(imax)
        return self
    
    ## Return the maximal number of CG iterations.
    # @return the maximal number of CG iterations.
    def getImax(self):
        return self.nMaxIterations
        
    ## Rises LinearSolverEvents.STARTING event
    def starting(self):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.STARTING)
    
    ## Rises LinearSolverEvents.CALC_STARTING event
    def calcStarting(self, ):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.CALC_STARTING)
    
    ## Rises LinearSolverEvents.ITERATION_COMPLETE event
    def iterationComplete(self, ):
        LinearSolver.notifyEventControllers(self, LinearSolverEvents.ITERATION_COMPLETE)
    
    ## Rises LinearSolverEvents.COMPLETE event
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
                serializationString += "'" + attrName + "'" + " : '" + attrValue + "',\n"
                
        serializationString = "{" + serializationString.rstrip(",\n") + "}\n"
        return serializationString    
    
    
    ## Restores the CGSolver object from the json object with attributes.
    #
    # @param cls python keyword (do not specify)
    # @param jsonObject A json object.
    # @return The restored SGSolver object.
    @classmethod
    def fromJson(cls, jsonObject):
        cg = CGSolver()
#        if jsonObject.has_key('accuracy'):
#            cg.setEpsilon( jsonObject['accuracy'] )
#        if jsonObject.has_key('imax'):
#            cg.setImax( jsonObject['imax'] )
        if jsonObject.has_key('delta_0'):
            cg.delta_0 = jsonObject['delta_0']
        if jsonObject.has_key('delta_new'):
            cg.delta_new = jsonObject['delta_new']
        if jsonObject.has_key('nIterations'):
            cg.nIterations = jsonObject['nIterations']
        if jsonObject.has_key('nMaxIterations'):
            cg.nMaxIterations = jsonObject['nMaxIterations']
        if jsonObject.has_key('myEpsilon'):
            cg.setEpsilon( jsonObject['myEpsilon'] )
        if jsonObject.has_key('residuum'):
            cg.residuum = float(jsonObject['residuum'])
        return cg



# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

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


## Abstract class of Subscribers of LinearSolverEvents. The classes that wants to obtain
# the progress notifications from LinearSolver subclasses should implement this class. See @link
# bin.learner.solver.LinearSolver.LinearSolver documentation of Learner@endlink for details.
class SolverEventController(object):

    ##
    #Handles events from LinearSolver 
    #
    #@param subject: LinearSolver object
    #@param event: Event Status of type LinearSolverEvents
    ##  
    def handleSolvingEvent(self, subject, event):
        raise NotImplementedError()
    
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.
    def toString(self):
        return "{'module' : '" + self.__module__ + "'}\n"
    
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.
    # @todo (khakhutv) write a test for this method
    def __repr__(self):
        return '{' + self.toString().lstrip("{").rstrip("}\n") + '}'

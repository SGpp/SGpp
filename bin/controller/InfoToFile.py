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

from bin.controller.InfoToScreen import InfoToScreen
import sys


## The class processes the progress information from Learner and LinearSolver and
# stores it into the file.
class InfoToFile(InfoToScreen):
    
    ## Filename, where output should be written
    filename = None 
    
    ##
    #Constructor
    #
    #@param filename: filename where output should be written as string 
    def __init__(self, filename):
        self.filename = filename
        
        
    ##
    #Handles events from Linear Solver 
    #
    #@param subject: Linear Solver object
    #@param status: Event Status of type LinearSolverEvents
    ##       
    def handleSolvingEvent(self, subject, status):
        file = open(self.filename, "a")
        tmpout = sys.stdout
        sys.stdout = file
        InfoToScreen.handleSolvingEvent(self, subject, status)
        sys.stdout = tmpout
        file.close()
        
        
    ##
    #Handles events from Learner 
    #
    #@param subject: Learner object
    #@param status: Event Status of type LearnerEvents
    ##        
    def handleLearningEvent(self, subject, status):
        file = open(self.filename, "a")
        tmpout = sys.stdout
        sys.stdout = file
        InfoToScreen.handleLearningEvent(self, subject, status)
        sys.stdout = tmpout
        file.close()
        
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.  
    def toString(self):
        serializationString = InfoToScreen.toString(self).rstrip().lstrip("{").rstrip("}")
        serializationString += ", 'filename' : '" + self.filename + "'"
        serializationString = "{" + serializationString + "}\n"
        return serializationString
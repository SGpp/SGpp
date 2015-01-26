// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

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

from pysgpp import DataVector


## Class keeps all information, which was learned during the learning process.
# Currently, only the alpha vector is stored, by in the future more information may come.
class LearnedKnowledge(object):

    
    __alphas = None         #DataVector with learned alpha
    
    
    ## Constructor
    def __init__(self):
        self.__alphas = DataVector(1) #a dummy DataVector object
        

    ##Restores the state which is saved in the given memento
    #
    #@param memento the memento object
    def setMemento(self, memento):
        self.update(memento)
    
    
    ##Creates a new memento to hold the current state
    #
    #@return a new memento
    def createMemento(self):
        # @todo (khakhutv) maybe a better Memento as just DataVector object required?
        return self.__alphas


    ##Returns current alpha vector
    #
    #@return: DataVector of current alpha
    def getAlphas(self):
        return self.__alphas


    ##Alters the current alpha vector
    #
    #@param alpha: new alpha vector
    #@return: LearnedKnowledge itself
    def update(self, alpha):
        self.__alphas = alpha
        return self





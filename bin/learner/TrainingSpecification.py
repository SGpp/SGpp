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

## @package TrainingSpecification
# @ingroup learner
# @brief Stores different parameters of the learning process
# @version $CURR$

class TrainingSpecification(object):
    """ generated source for TrainingSpecification

    """
    __adaptPoints = 0
    __l = None
    __adaptRate = 0
    __cOperator = None
    __bOperator = None

    def setAdaptPoints(self, value):
        self.__adaptPoints = value


    def setL(self, value):
        self.__l = value


    def setAdaptRate(self, value):
        self.__adaptRate = value


    def setCOperator(self, value):
        self.__cOperator = value


    def setBOperator(self, value):
        self.__bOperator = value


    def getAdaptPoints(self):
        return self.__adaptPoints


    def getL(self):
        return self.__l


    def getAdaptRate(self):
        return self.__adaptRate


    def getCOperator(self):
        return self.__cOperator


    def getBOperator(self):
        return self.__bOperator
    
    def getNumOfPointsToRefine(self, refinablePoints):
        ratePoints = self.__adaptRate * refinablePoints
        if self.__adaptPoints == 0:
            return ratePoints
        elif self.__adaptRate == 0:
            return self.__adaptPoints
        else: 
            return min(ratePoints, self.__adaptPoints)





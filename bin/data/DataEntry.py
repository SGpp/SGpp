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
from DataSpecification import DataSpecification


## A container for tuple of a point and its corresponding value
class DataEntry(object):

    ## DataVector for value
    value = None        
    
    ## DataVector for point
    point = None        
    
    
    ##Constructor
    #
    # @param point: DataVector data point
    # @param value: DataVector value of function in the point
    def __init__(self,point, value):
        self.point = point
        self.value = value
    
    
    ## Returns the value
    #
    # @return: DataVector value
    def getValue(self):
        return self.value
    
    
    ## Returns the data point
    #
    # @return: DataVector point
    def getPoint(self):
        return self.point

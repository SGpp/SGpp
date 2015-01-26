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


## Abstract class defines the interface for storing and loading of input data.
class DataAdapter(object):


    ## Store data into file
    # as it is an abstract class, this function is not implemented!
    # 
    # @param points: DataVector with points
    # @param values: DataVector with values, default None
    # @param attributes: dictionary with attributes of dataset, default None
    def save(self, points, values = None, attributes = None):
        raise NotImplementedError


    ## Reads dataset from file
    # as it is an abstract class, this function is not implemented!
    #
    # @param name: String for category of data set (train or test), default "train"
    # @return DataContainer with data set
    def loadData(self, name="train"):
        raise NotImplementedError


    ## Loads attribute specification from file
    # as it is an abstract class, this function is not implemented!
    #
    # @return dictionary with attribute specification
    def loadSpecification(self,):
        raise NotImplementedError



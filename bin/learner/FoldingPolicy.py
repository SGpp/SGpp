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


import math

## Abstract class for providing functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs
# @todo (khakhutv) why is it class abstract?
class FoldingPolicy(object):


    level = None        #Folding level
    size = None         #Size of dataset
    dataset = None      #Dataset
    dataFold = []       #List of partitioned data sets
    window = None       #Number of points in one subset
    seq = None          #Sequence of indices of points from data set
    
    
    ##Constructor
    #
    #@param dataset: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    def __init__(self, dataset, level=1):
        self.dataFold = []
        self.level = level
        self.dataset = dataset
        self.size = dataset.getPoints().getSize()
        self.window = int( math.ceil( float(self.size) / self.level ) ) #number of points in validation set

        
    ##Implementation of iterator method next()
    #
    # @return: the next subset
    def next(self):
        for step in xrange(self.level):
            yield self.dataFold[step]
        return
    
    ## Create fold new data set
    # Brings points given by validationIndeces together as test subset and the rest of points
    # as train subset
    #
    # @param dataContainer: DataContainer with points
    # @param validationIndeces: list of indices for validation subset
    # @return: DataContainer partitioned data set
    def createFoldsets(self, dataContainer, validationIndeces):
        foldContainerValidation = dataContainer.getDataSubsetByIndexList(validationIndeces, "test")
        trainIndeces = [i for i in self.seq if i not in validationIndeces]
        foldContainerTrain = dataContainer.getDataSubsetByIndexList(trainIndeces, "train")
        return foldContainerTrain.combine(foldContainerValidation)

    
    
    ##Implementation of iterator method __iter__()
    # iterates through subsets
    def __iter__(self):
        return self.next()
            



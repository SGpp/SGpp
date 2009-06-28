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

## @package SequentialFoldingPolicy
# @ingroup learner
# @brief Folding rule
# @version $CURR$

import time
import random
import math

from FoldingPolicy import FoldingPolicy
from bin.data.ARFFAdapter import ARFFAdapter


class SequentialFoldingPolicy(FoldingPolicy):
    
    random = None
    seed = None
    seq = None
    window = None
    __dataFold = []
    
    def __init__(self, dataContainer, level):
        self.__dataFold = []
        FoldingPolicy.__init__(self,  dataContainer, level)
        self.window = int( math.ceil( float(self.size) / self.level ) ) #number of points in validation set
        self.seq = range(self.size)
        for step in xrange(self.level):
            validationIndeces = self.seq[ step * self.window : min((step+1) * self.window, self.size)]
            self.__dataFold.append(self.__createFoldsets(dataContainer, validationIndeces))
            

    def __createFoldsets(self, dataContainer, validationIndeces):
        foldContainerValidation = dataContainer.getDataSubsetByIndexList(validationIndeces, "test")
        trainIndeces = [i for i in self.seq if i not in validationIndeces]
        foldContainerTrain = dataContainer.getDataSubsetByIndexList(trainIndeces, "train")
        return foldContainerTrain.combine(foldContainerValidation)

            
        
    def next(self):
        for step in xrange(self.level):
            yield self.__dataFold[step]
        return
    
    def __iter__(self):
        return self.next()
            

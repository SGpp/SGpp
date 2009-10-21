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



import time
import random
import math

from FoldingPolicy import FoldingPolicy

## Provides functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs randomly
class RandomFoldingPolicy(FoldingPolicy):
    
    
    random = None   #Random number generator

    
    ##Constructor
    #
    #@param dataContainer: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    #@param seed: Integer seed, default None so it is set to the timestamp
    def __init__(self,  dataContainer, level=1, seed = None):
        FoldingPolicy.__init__(self,  dataContainer, level)
        self.window = int( math.ceil( self.size / self.level ) )
        if seed == None:
            self.seed = int(time.time())
        else:
            self.seed = seed
        self.random = random.seed(self.seed)
        self.seq = range(self.size)
        random.shuffle(self.seq, self.random)
        for step in xrange(self.level):
            validationIndeces = self.seq[ step * self.window : min((step+1) * self.window, self.size)]
            self.dataFold.append(self.createFoldsets(dataContainer, validationIndeces))


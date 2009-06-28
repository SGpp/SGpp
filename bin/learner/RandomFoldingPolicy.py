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

## @package RandomFoldingPolicy
# @ingroup bin.learner
# @brief Folding rule
# @version $CURR$

import time
import random
import math

from FoldingPolicy import FoldingPolicy


class RandomFoldingPolicy(FoldingPolicy):
    
    random = None   #Random number generator
    seed = None     #Seed, the random number generator is initialized with
    seq = None      #Sequence of indeces of points from data set
    window = None   #Number of points in one subset
    
    
    ##Constructor
    #
    #@param dataset: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    #@param seed: Integer seed, default None so it ist set to the timestamp
    def __init__(self,  dataset, level=1, seed = None):
        FoldingPolicy.__init__(self,  dataset, level)
        self.window = int( math.ceil( self.size / self.level ) )
        if seed == None:
            self.seed = int(time.time())
        else:
            self.seed = seed
        self.random = random.seed(self.seed)
        self.seq = range(self.size)
        random.shuffle(self.seq, self.random)


    ##Implementation of iterator method next()
    #
    # @return: the next subset     
    def next(self):
        for step in xrange(self.level):
            if step != (self.level-1):
                yield self.seq[ step * self.window : (step+1) * self.window ]
            else:
                yield self.seq[ step * self.window : ]
        return


    ##Implementation of iterator method __iter__()
    # iterates through subsets   
    def __iter__(self):
        return self.next()
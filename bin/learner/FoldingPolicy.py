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

## @package FoldingPolicy
# @ingroup bin.learner
# @brief Abstract class for folding rules
# @version $CURR$

class FoldingPolicy(object):


    level = None    #Folding level
    size = None     #Size of dataset
    dataset = None  #Dataset
    
    
    ##Constructor
    #
    #@param dataset: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    def __init__(self, dataset, level=1):
        self.level = level
        self.dataset = dataset
        self.size = dataset.getPoints().getSize()

        
    ##Implementation of iterator method next()
    #
    # @return: the next subset
    def next(self):
        for step in xrange(self.level):
                yield range(self.size)
        return
    
    
    ##Implementation of iterator method __iter__()
    # iterates through subsets
    def __iter__(self):
        return self.next()
            



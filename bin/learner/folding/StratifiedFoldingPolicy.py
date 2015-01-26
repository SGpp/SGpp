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

from FoldingPolicy import FoldingPolicy
import math

## Provides functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs with equal distribution
# of points from two different classes between folds.
# This class corresponds to the old doFoldr() and doFoldStratified() methods.
class StratifiedFoldingPolicy(FoldingPolicy):
    
    
    ##Constructor
    #
    #@param dataContainer: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    def __init__(self, dataContainer, level=1):
        FoldingPolicy.__init__(self,  dataContainer, level)
        values = dataContainer.getValues()
        self.seq = range(self.size)
        indecesOfPositives = []
        indecesOfNegatives = []
        for i in xrange(self.size):
            if values[i] < 0: indecesOfNegatives.append(i)
            else: indecesOfPositives.append(i)
        windowPos = int(math.floor(len(indecesOfPositives) / self.level))
        windowNeg = int(math.floor(len(indecesOfNegatives) / self.level))
        for step in xrange(self.level):
            validationIndeces = indecesOfPositives[ step * windowPos : ((step+1) * windowPos if (step+1) < self.level else len(indecesOfPositives))] + \
                                indecesOfNegatives[ step * windowNeg : ((step+1) * windowNeg if (step+1) < self.level else len(indecesOfNegatives))]
            self.dataFold.append(self.createFoldsets(dataContainer, validationIndeces))
            
            
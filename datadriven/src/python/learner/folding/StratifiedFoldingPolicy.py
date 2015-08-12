# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##############################################################################
                                    #
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
            
            
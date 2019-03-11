# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


import time
import math

from pysgpp.extensions.datadriven.learner.folding.FoldingPolicy import FoldingPolicy
from pysgpp.extensions.datadriven.data.ARFFAdapter import ARFFAdapter

## Provides functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs sequentially
# This class corresponds to the old doFolds() method.
class SequentialFoldingPolicy(FoldingPolicy):


    ##Constructor
    #
    #@param dataContainer: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    def __init__(self, dataContainer, level):
        FoldingPolicy.__init__(self,  dataContainer, level)
        self.seq = list(range(self.size))
        for step in range(self.level):
            validationIndeces = self.seq[ step * self.window : min((step+1) * self.window, self.size)]
            self.dataFold.append(self.createFoldsets(dataContainer, validationIndeces))

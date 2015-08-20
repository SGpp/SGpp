# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import time
import random
import math

from FoldingPolicy import FoldingPolicy

## Provides functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs randomly
# This class corresponds to the old doFold() method.
class RandomFoldingPolicy(FoldingPolicy):
    
    
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
        ## Random number generator
        self.random = random.seed(self.seed)
        self.seq = range(self.size)
        random.shuffle(self.seq, self.random)
        for step in xrange(self.level):
            validationIndeces = self.seq[ step * self.window : min((step+1) * self.window, self.size)]
            self.dataFold.append(self.createFoldsets(dataContainer, validationIndeces))

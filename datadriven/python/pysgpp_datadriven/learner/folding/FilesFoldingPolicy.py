# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

##############################################################################
                                    #
#############################################################################

from FoldingPolicy import FoldingPolicy
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
import math

## Provides functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs from the set of files
# This class corresponds to the old doFoldf() method.
class FilesFoldingPolicy(FoldingPolicy):


    ##Constructor
    #
    #@param dataContainer: DataContainer with data set
    #@param level: Integer folding level, default value: 1. This parameter is used for compatibility only. The folding level will be set to the number of files.
    def __init__(self, dataContainer, level=1):
        FoldingPolicy.__init__(self,  dataContainer, level)
        datasets = []
        fileCounter = 0
        # It is expected, that several files are stored in the data container
        # with category name "train0", "train1" etc. If the category name "trainN"
        # doesn't exist and getDataSubsetByCategory() raises an Exception, it
        # means, we've gathered all data sets.
        while True:
            try:
                datasets.append(dataContainer.getDataSubsetByCategory(DataContainer.TRAIN_CATEGORY + str(fileCounter)))
                fileCounter += 1
            except:
                break

        # as in the old doFoldf() method the folding level is determined by the
        # number of files:
        self.level = fileCounter

        for step in xrange(self.level):
            self.dataFold.append(
                                 # merge all data containers except the one with number =step
                                 DataContainer.merge(datasets[:step] + datasets[step+1:])
                                 )

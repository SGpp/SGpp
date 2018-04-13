# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import math


## Abstract class for providing functionality for accomplishment of learning with cross-validation
# by generating a set of training data/validation data pairs
class FoldingPolicy(object):
    
    ##Constructor
    #
    #@param dataset: DataContainer with data set
    #@param level: Integer folding level, default value: 1
    def __init__(self, dataset, level=1):
        ##List of partitioned data sets
        self.dataFold = []
        
        ##Folding level
        self.level = level
        
        ##Dataset
        self.dataset = dataset
        
        ##Size of dataset
        self.size = dataset.getPoints().getNrows()
        
        ##Number of points in one subset
        self.window = int( math.ceil( float(self.size) / self.level ) ) #number of points in validation set
        
        ##Sequence of indices of points from data set
        self.seq = None 

        
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
            


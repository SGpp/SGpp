# This file is part of sgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#

from pysgpp import *

class DataProvider(object):
    def __init__(self, data, mode):
        self.__data = data
        self.__nextMethod = self.data_providers[mode]["split"]
        self.__constructMethod = self.data_providers[mode]["construct"]
    
    def __iter__(self):
        return self.__nextMethod(self, self.data)
    
    def construct(self):
        return self.__constructMethod(self, self.data)

    def getData(self):
        return self.__data

    
    def constructNormal(self, data):
        if len(data) != 1:
            raise Exception("Only one data file supported.")
    
        return (data[0], len(data[0]["data"]))
    
    def splitNormal(self, data):
        yield (self.buildTrainingVector(data), self.buildYVector(data)), None
        return
    
    ## Random Fold
    def constructFold(self, data):
        if len(data) != 1:
            raise Exception("Only one data file supported.")
        
        return split_n_folds(data[0], options.f_level, options.seed), len(data[0]["data"])
    
    ## Stratified Fold
    def constructFoldr(self, data):
        if len(data) != 1:
            raise Exception("Only one data file supported.")
        
        return split_n_folds_stratified(data[0], options.f_level, options.seed), len(data[0]["data"])
    
    
    ## Sequential Fold
    def constructFolds(self, data):
        if len(data) != 1:
            raise Exception("Only one data file supported.")
        
        return split_n_folds_sequential(data[0], options.f_level), len(data[0]["data"])
    
    
    def splitFold(self, data):
        dvec, cvec = data
        for i in xrange(options.f_level):
            training = assembleTrainingVector(dvec,cvec,i)
            testing = (dvec[i], cvec[i])
            yield (training, testing)
        return
    
    def buildTrainingVector(self, data):
        dim = len(data["data"])
        training = DataVector(len(data["data"][0]), dim)
        
        for i in xrange(len(data["data"][0])):
            for d in xrange(dim):
                training[i*dim + d] = data["data"][d][i]
        
        return training
    
    def buildYVector(self, data):
        y = DataVector(len(data["classes"]))
        for i in xrange(len(data["classes"])):
            y[i] = data["classes"][i]
            
        return y
    
    #-------------------------------------------------------------------------------
    # For n-fold-cv:
    # assemble training vector for fold <omit>
    #-------------------------------------------------------------------------------
    def assembleTrainingVector(self, dvecs, cvecs, omit):
        size = 0
        for dataset in dvecs:
            size = size + dataset.getSize()
        
        size = size - dvecs[omit].getSize()
        
        training = DataVector(size, dvecs[0].getDim())
        classes = DataVector(size)
        
        i=0
        for vec in dvecs:
            if vec == dvecs[omit]:
                continue
            for x in xrange(len(vec)):
                training[i] = vec[x]
                i = i + 1
        
        i=0
        for vec in cvecs:
            if vec == cvecs[omit]:
                continue 
            for x in xrange(len(vec)):
                classes[i] = vec[x]
                i = i + 1
                
        return training, classes
    
    ## list of data providers
    # construct should bring the data in a suitable form.
    # split should be a generator returning a tuple (training, testing) for the current fold.
    # training/testing should be of the form (data, classes) or None if not present
    data_providers = {
            "normal" : {"construct" : constructNormal, "split" : splitNormal},
            "fold" : {"construct" : constructFold, "split" : splitFold},
            "folds" : {"construct" : constructFolds, "split" : splitFold},                  
            "foldr" : {"construct" : constructFoldr, "split" : splitFold},                  
            }

    data = property(getData, "Data's Docstring")

    

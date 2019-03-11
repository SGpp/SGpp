# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

## Abstract class defines the interface for storing and loading of input data.
class DataAdapter(object):


    ## Store data into file
    # as it is an abstract class, this function is not implemented!
    # 
    # @param points: DataVector with points
    # @param values: DataVector with values, default None
    # @param attributes: dictionary with attributes of dataset, default None
    def save(self, points, values = None, attributes = None):
        raise NotImplementedError


    ## Reads dataset from file
    # as it is an abstract class, this function is not implemented!
    #
    # @param name: String for category of data set (train or test), default "train"
    # @return DataContainer with data set
    def loadData(self, name="train"):
        raise NotImplementedError


    ## Loads attribute specification from file
    # as it is an abstract class, this function is not implemented!
    #
    # @return dictionary with attribute specification
    def loadSpecification(self,):
        raise NotImplementedError


# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

from pysgpp import DataVector
from DataSpecification import DataSpecification

## A container for tuple of a point and its corresponding value
class DataEntry(object):

    ## DataVector for value
    value = None        
    
    ## DataVector for point
    point = None        
    
    
    ##Constructor
    #
    # @param point: DataVector data point
    # @param value: DataVector value of function in the point
    def __init__(self,point, value):
        self.point = point
        self.value = value
    
    
    ## Returns the value
    #
    # @return: DataVector value
    def getValue(self):
        return self.value
    
    
    ## Returns the data point
    #
    # @return: DataVector point
    def getPoint(self):
        return self.point
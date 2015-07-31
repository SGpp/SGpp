# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

from time import time
import random

## The class encapsulates the description of the data set, like names of dimensions
# and source file.
class DataSpecification(object):

    
    ## Checks if the data set was already saved as a file
    # @return: True, if data set was saved, otherwise: False
    def isSaved(self):
        return self.__isSaved

    
    ## Marks the data set as saved
    def setSaved(self):
        self.__isSaved = True

    ## The file name
    filename = None
    
    ## Dictionary of attributes
    attributes = None
    
    __isSaved = None
    
    
    ## Constructor
    def __init__(self):
        self.filename = ''
        self.attributes = {}
        self.__isSaved = False
    
    
    ## Sets the file name
    # @param filename: new file name
    def setFilename(self, filename):
        self.filename = filename
    
    
    ## Returns the file name
    # @return: the file name
    def getFilename(self):
        #if filename is not set, it means the data subset was generated dynamically
        #so create new name, don't forget to save the data into the corresponding file!
        if self.filename == '':
            self.setFilename(self.generateFilename())
        return self.filename
    
    
    ## Generates a new random file name
    # @return: new random file name
    def generateFilename(self):
        rnd = str(int(random.random()*1000)) #random post-fix
        return 'data.' + str(time()) + rnd + '.arff'
    
    ## Adds an attribute to the attribute collection
    # @param attributeName: attribute name
    # @param attributeValue: attribute value
    def addAttribute(self, attributeName, attributeValue):
        self.attributes[attributeName] = attributeValue
    
    
    ## Returns the collection of attributes stored in the specification
    # @return: collection of attributes    
    def getAttributes(self):
        return self.attributes
    
    
    ## Deletes all data from attributes collection
    def cleanAttributes(self):
        self.attributes = {}
    
    ## Fills data specification with a number of numeric attributes with names x0, x1, x2...  
    # @param dim the number of attributes
    def createNumericAttributes(self, dim):
        self.cleanAttributes()
        for i in xrange(0,dim):
            self.addAttribute('x'+str(i), 'NUMERIC')
        self.addAttribute('class', 'NUMERIC')
    
    
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.
    def toString(self):
        serializationString = "'module' : '" + self.__module__ + "',\n"
        serializationString += '"filename":"' + self.getFilename() + '"\n'
        for key, value in self.getAttributes().items():
            serializationString = serializationString + ',"' + key + '" : "' + value + '"'
        return '{' + serializationString + '}'


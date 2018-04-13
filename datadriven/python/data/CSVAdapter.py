# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#############################################################################
                                    #
#############################################################################
from DataSpecification import DataSpecification


import re
import gzip
import csv
from DataAdapter import DataAdapter
from pysgpp import DataVector, DataMatrix

from DataContainer import DataContainer


## Class implements the interface of DataAdapter for storing and restoring of input
# data into / from files in CSV-format.
class CSVAdapter(DataAdapter):

    ## Filename associated with data
    filename = None


    ## Constructor
    #
    # @param filename: Filename as String
    def __init__(self, filename=""):
        self.filename = filename
    
    
    ## Store data into file
    # 
    # @param points: DataVector with points
    # @param values: DataVector with values, default None
    # @param attributes: dictionary with attributes of dataset, default None
    def save(self, points, values = None, attributes = None):

        fout = self.__gzOpen(self.filename, "w")
        dim = points.getNcols()
        size = points.getNrows()
        point = DataVector(dim)
        
        fout.write("@RELATION \"%s\"\n\n" % self.filename)
        
        hasclass = False
        if values != None:
            hasclass = True
        
        if attributes == None:
            for i in xrange(dim):
                fout.write("@ATTRIBUTE x%d NUMERIC\n" % i)
               
            if hasclass:
                fout.write("@ATTRIBUTE class NUMERIC\n")
        else:
            for key in attributes.keys():
                fout.write("@ATTRIBUTE %s %s\n" % (key, attributes[key]))
            
        fout.write("\n@DATA\n")
            
        if attributes != None:
            fout.write("%s\n"%','.join(attributes))

        for row in xrange(size):
            points.getRow(row, point)
            lout = []
            for col in xrange(dim):
                lout.append(point[col])
            if hasclass:
                lout.append(values[row])
            string = ",".join(str(i) for i in lout)
            fout.write(string+"\n")

        fout.close()


    ## Reads dataset from file
    #
    # @param name: String for category of data set (train or test), default "train"
    # @param delimiter optional delimiter character. Default: ','
    # @param target_col optional number of target column. Default: -1
    # @return DataContainer with data set
    def loadData(self, name = "train", delimiter=',', target_col=-1):
        fin = self.__gzOpen(self.filename, "r")
        reader = csv.reader(fin, delimiter = delimiter)
        data = []
        classes = []
        hasclass = False
        target_col = -1
        
        first_line = reader.next()
        
        # training set has to contain targets
        if name == 'train':
            if len(first_line) <= target_col: 
                raise Exception('Target column does not match total column number.')
            for i in xrange(len(first_line)-1):
                data.append([])
            hasclass = True
        
        # test set may contain target values
        if name == 'test':
            if len(first_line) > target_col: 
                for i in xrange(len(first_line)-1):  data.append([])
                hasclass = True
            else:
                for i in xrange(len(first_line)):  data.append([])
                hasclass = False
        
        
        # skip header if available
        if first_line[0].isalpha():
            pass
        else:
            line = first_line
            if hasclass:
                classes.append(float(line[target_col]))
                line.remove(line[target_col])
            for i in xrange(len(line)):
                data[i].append(float(line[i]))
        
        for line in reader:
            if hasclass:
                classes.append(float(line[target_col]))
                line.remove(line[target_col])
            for i in xrange(len(line)):
                data[i].append(float(line[i]))
                
        # cleaning up and return
        fin.close()
        
        dim = len(data)
        size = len(data[0])
        dataMatrix = DataMatrix(size, dim)
        tempVector = DataVector(dim)
        valuesVector = DataVector(size)
        for rowIndex in xrange(size):
            for colIndex in xrange(dim):
                tempVector[colIndex] = data[colIndex][rowIndex]
            dataMatrix.setRow(rowIndex, tempVector)
            valuesVector[rowIndex] = classes[rowIndex]
            
        return DataContainer(dataMatrix, valuesVector, name, self.filename)


    ## Loads attribute specification from file
    #
    # @return dictionary with attribute specification
    def loadSpecification(self,  delimiter=','):
        spec = DataSpecification()
        spec.setFilename(self.filename)
        fin = self.__gzOpen(self.filename, "r")
        reader = csv.reader(fin, delimiter = delimiter)
        first_line = reader.next()
        second_line = reader.next()
        if first_line[0].isalpha():
            for col in xrange(len(first_line)):
                try:
                    i = float(second_line[col])
                except ValueError, TypeError:
                    # not numeric
                    spec.addAttribute(first_line[col], "string")
                else:
                    # numeric
                    spec.addAttribute(first_line[col], "numeric")
        return spec
    

    ## Opens a file. If the file ends with ".gz", automatically gzip compression
    # is used for the file. Returns the filedescriptor
    # @param filename
    # @param mode default: "r" for read only
    # @return file descriptor
    def __gzOpen(self, filename, mode="r"):
        # gzip-file
        if re.match(".*\.gz$", filename):
            # mode set for binary data?
            if not mode[-1] == "b":
                mode += "b"
            fd = gzip.open(filename, mode)
        # non gzip-file
        else:
            fd = open(filename, mode)
        return fd


# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

#############################################################################
                                    #
#############################################################################
from pysgpp.extensions.datadriven.data.DataSpecification import DataSpecification


import re
import gzip
from pysgpp.extensions.datadriven.data.DataAdapter import DataAdapter
from pysgpp import DataVector, DataMatrix

from pysgpp.extensions.datadriven.data.DataContainer import DataContainer

## Class implements the interface of DataAdapter for storing and restoring of input
# data into / from files in <a href="http://www.cs.waikato.ac.nz/~ml/weka/arff.html" 
# target="new">ARFF-format</a>.
class ARFFAdapter(DataAdapter):

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
        
        fout.write(b"@RELATION \"%s\"\n\n" % self.filename.encode())
        
        hasclass = False
        if values != None:
            hasclass = True
        
        if attributes == None:
            for i in range(dim):
                fout.write(b"@ATTRIBUTE x%d NUMERIC\n" % i)
               
            if hasclass:
                fout.write(b"@ATTRIBUTE class NUMERIC\n")
        else:
            for key in list(attributes.keys()):
                fout.write(b"@ATTRIBUTE %s %s\n" % (key.encode(), attributes[key].encode()))
            
        fout.write(b"\n@DATA\n")
            

        for row in range(size):
            points.getRow(row, point)
            lout = []
            for col in range(dim):
                lout.append(point[col])
            if hasclass:
                lout.append(values[row])
            string = ",".join(str(i) for i in lout)
            fout.write((string+"\n").encode())

        fout.close()


    ## Reads dataset from file
    #
    # @param name: String for category of data set (train or test), default "train"
    # @return DataContainer with data set
    def loadData(self, name = "train"):
        fin = self.__gzOpen(self.filename, "r")
        data = []
        classes = []
        hasclass = False
    
        # get the different section of ARFF-File
        for line in fin:
            sline = line.strip().lower()
            if sline.startswith(b"%") or len(sline) == 0:
                continue
    
            if sline.startswith(b"@data"):
                break
            
            if sline.startswith(b"@attribute"):
                value = sline.split()
                if value[1].startswith(b"class"):
                    hasclass = True
                else:
                    data.append([])
        
        #read in the data stored in the ARFF file
        for line in fin:
            sline = line.strip()
            if sline.startswith(b"%") or len(sline) == 0:
                continue
    
            values = sline.split(b",")
            if hasclass:
                classes.append(float(values[-1]))
                values = values[:-1]
            for i in range(len(values)):
                data[i].append(float(values[i]))
                
        # cleaning up and return
        fin.close()
        
        dim = len(data)
        size = len(data[0])
        dataMatrix = DataMatrix(size, dim)
        tempVector = DataVector(dim)
        valuesVector = DataVector(size)
        for rowIndex in range(size):
            for colIndex in range(dim):
                tempVector[colIndex] = data[colIndex][rowIndex]
            dataMatrix.setRow(rowIndex, tempVector)
            valuesVector[rowIndex] = classes[rowIndex]
            
        return DataContainer(points=dataMatrix, values=valuesVector, name=name, filename=self.filename)


    ## Loads attribute specification from file
    #
    # @return dictionary with attribute specification
    def loadSpecification(self):
        spec = DataSpecification()
        spec.setFilename(self.filename)
        fin = self.__gzOpen(self.filename, "r")
        for line in fin:
            sline = line.strip().lower()
            # comments:
            if sline.startswith(b"%") or len(sline) == 0:
                continue
            
            # attributes:
            elif sline.startswith(b"@attribute"):
                attrLine = line.strip().split()
                spec.addAttribute(attrLine[1], attrLine[2])
            
            # data block therefore no attributes will come
            elif sline.startswith(b"@data"):
                break
            
        return spec
    

    ## Opens a file. If the file ends with ".gz", automatically gzip compression
    # is used for the file. Returns the filedescriptor
    # @param filename
    # @param mode default: "r" for read only
    # @return file descriptor
    def __gzOpen(self, filename, mode="r"):
        # mode set for binary data?
        if not mode[-1] == "b":
            mode += "b"
        # gzip-file
        if re.match(".*\.gz$", filename):
            fd = gzip.open(filename, mode)
        # non gzip-file
        else:
            fd = open(filename, mode)
        return fd


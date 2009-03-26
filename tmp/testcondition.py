## This is alex's test file

from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from painlesscg import cg,sd,cg_new
from math import sqrt
import random
import numpy as np

from array import array

try:
    import psyco
    psyco.full()
    print "Using psyco"
except:
    pass


#-------------------------------------------------------------------------------
## Outputs a deprecated warning for an option
# @param option Parameter set by the OptionParser
# @param opt Parameter set by the OptionParser
# @param value Parameter set by the OptionParser
# @param parser Parameter set by the OptionParser
def callback_deprecated(option, opt, value, parser):
    print "Warning: Option %s is deprecated." % (option)
    

## Reads in an ARFF file:
# The data is stored in lists. There is a value list for every dimension of the data set. e.g. 
# [[2, 3],[1, 1]] are the data points P_1(2,1) and P_2(3,1)
#
# @param filename the file's filename that should be read
# @return returns a set of a array with the data (named data), a array with the classes (named classes) and the filename named as filename
def readDataARFF(filename):
    fin = open(filename, "r")
    data = []
    classes = []
    hasclass = False

    # get the different section of ARFF-File
    for line in fin:
        sline = line.strip().lower()
        if sline.startswith("%") or len(sline) == 0:
            continue

        if sline.startswith("@data"):
            break
        
        if sline.startswith("@attribute"):
            value = sline.split()
            if value[1].startswith("class"):
                hasclass = True
            else:
                data.append([])
    
    #read in the data stored in the ARFF file
    for line in fin:
        sline = line.strip()
        if sline.startswith("%") or len(sline) == 0:
            continue

        values = sline.split(",")
        if hasclass:
            classes.append(float(values[-1]))
            values = values[:-1]
        for i in xrange(len(values)):
            data[i].append(float(values[i]))
            
    # cleaning up and return
    fin.close()
    return {"data":data, "classes":classes, "filename":filename}

    
#-------------------------------------------------------------------------------
## Builds the training data vector
# 
# @param data a list of lists that contains the points a the training data set, coordinate-wise
# @return a instance of a DataVector that stores the training data
def buildTrainingVector(data):
    dim = len(data["data"])
    training = DataVector(len(data["data"][0]), dim)
    
    # i iterates over the data points, d over the dimension of one data point
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            training[i*dim + d] = data["data"][d][i]
    
    return training


def openFile(filename):
    try:
        data = readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        sys.exit(1)
        
    if data.has_key("classes") == False:
        print ("No classes found in the given File " + filename + "!")
        sys.exit(1)
        
    return data


def generateCMatrix(factory, level, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
    
    laplace = factory.createOperationLaplace()
    
    # create vector
    alpha = DataVector(storage.size())
    erg = DataVector(storage.size())
    col = 0

    # create stiffness matrix
    m = np.zeros( (storage.size(), storage.size()) )

    for i in xrange(storage.size()):
        # apply unit vectors
        alpha.setAll(0)
        alpha[i] = 1
        laplace.mult(alpha, erg)
        if verbose:
            print erg, erg.sum()
        
        #Sets the column in m
        for j in xrange(storage.size()):
            m[j,col] = erg[j]
            
        col = col + 1

    return m


def generateBMatrix(factory, level, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    # create vector
    alpha = buildTrainingVector(openFile("function.out"))
    point = DataVector(1, storage.size())
    erg = DataVector(1, storage.dim())
    erg.setAll(0.0)
    col = 0
    # create B matrix
    m = np.zeros( (storage.dim(), storage.size()) )

    for i in xrange(storage.size()):
        # apply unit vectors
        point.setAll(0)
        point[i] = 1
        b.multTranspose(alpha, point, erg)
        if verbose:
            print erg, erg.sum()

        #Sets the column in m
        for j in xrange(storage.dim()):
            m[j,col] = erg[j]

        col = col + 1
        
    return m


def MapleMatrixString(m, Name, r, c):
    maplematrix = Name + ":=matrix(" + str(r) + "," + str(c) + ", ["
    
    print r
    print c

    for i in range(r):
        for j in range(c):
            if j < (c-1) or i < (r-1):
                maplematrix = maplematrix + str(round(m[i*r + j],10)) + ","
            else:
                maplematrix = maplematrix + str(round(m[i*r + j],10))
    
    maplematrix = maplematrix + "]):"
        
    return maplematrix


def build_DM_Matrices():
    factory = Grid.createLinearGrid(2)
    level = 7
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    aem = 194
    lam = 0.01
    
    print "generating laplacian matrix..."
    laplace_m = generateCMatrix(factory, level)
    print laplace_m
    print "generating transposed B matrix..."
    Bt_m = generateBMatrix(factory, level)
    print Bt_m
    print "generating B from transposed B..."
    B_m = Bt_m.transpose()
    print B_m
    print "multiplying B and B^T..." 
    B_res = np.dot(B_m,Bt_m)
    print B_res
    print "multiplying aem*lambda*C..."
    C = aem * lam * laplace_m
    print "adding C and B_res..."
    C = C + B_res
    print C
    print "calculating condition number..."
    cond = np.linalg.cond(C)
    
    print cond
    
    #write maple file
    #print "started writing maple file..."
    #fout = file("testcondition.maple", "w")
    #fout.write("with(linalg):\n")
    #fout.write(MapleMatrixString(laplace_m, "C", factory.getStorage().size(), factory.getStorage().size()))
    #fout.write("\n")
    #fout.write(MapleMatrixString(B_m, "B", factory.getStorage().size(), factory.getStorage().dim()))
    #fout.close()
    print "done!"
    
    
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    build_DM_Matrices()

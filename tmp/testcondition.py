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
        
        #Sets the column in m
        for j in xrange(storage.size()):
            m[j,col] = erg[j]
            
        col = col + 1

    return m


def generateBBTMatrix(factory, level, training, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    alpha = DataVector(storage.size())
    temp = DataVector(training.getSize())
    
    erg = DataVector(alpha.getSize())
    
    col = 0
    
    # create B matrix
    m = np.zeros( (storage.size(), storage.size()) )
    
    #print training

    for i in xrange(storage.size()):
        # apply unit vectors
        temp.setAll(0.0)
        erg.setAll(0.0)
        alpha.setAll(0.0)
        alpha[i] = 1.0
        b.multTranspose(alpha, training, temp)
        b.mult(temp, training, erg)
        
        #Sets the column in m
        for j in xrange(storage.size()):
            m[j,col] = erg[j]

        col = col + 1
        
    return m


#def MapleMatrixString(m, Name, r, c):
#    maplematrix = Name + ":=matrix(" + str(r) + "," + str(c) + ", ["
#    
#    print r
#    print c

#    for i in range(r):
#        for j in range(c):
#            if j < (c-1) or i < (r-1):
#                maplematrix = maplematrix + str(round(m[i*r + j],10)) + ","
#            else:
#                maplematrix = maplematrix + str(round(m[i*r + j],10))
    
#    maplematrix = maplematrix + "]):"
        
#    return maplematrix


def build_DM_Matrices():
    factory = Grid.createLinearBoundaryUScaledGrid(2)
    level = 3
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    training = buildTrainingVector(openFile('twospirals.wieland.arff.gz'))
    #training = buildTrainingVector(openFile('data_dim_1_nops_8_float.arff.gz'))
    
    aem = 194
    lam = 0.0001
    
    print "generating laplacian matrix..."
    laplace_m = generateCMatrix(factory, level)
    print laplace_m
    print "generating B*B^T matrix..."
    B_res = generateBBTMatrix(factory, level, training) #np.dot(B_m,Bt_m)
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
    #print "done!"
    
    
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    build_DM_Matrices()

## This is alex's test file

from optparse import OptionParser
import sys
sys.path.append("../bin")


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


def generateBBTMatrix(factory, training, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    alpha = DataVector(storage.size())
    temp = DataVector(training.getSize())
    
    erg = DataVector(alpha.getSize())
    
    col = 0
    
    # create B matrix
    m = np.zeros( (storage.size(), storage.size()) )
    
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


def generateBTMatrix(factory, training, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    alpha = DataVector(storage.size())
    temp = DataVector(training.getSize())
        
    col = 0
    
    # create BT matrix
    m = np.zeros( (training.getSize(), storage.size()) )
    
    for i in xrange(storage.size()):
        # apply unit vectors
        temp.setAll(0.0)
        alpha.setAll(0.0)
        alpha[i] = 1.0
        b.multTranspose(alpha, training, temp)
        
        #Sets the column in m
        for j in xrange(training.getSize()):
            m[j,col] = temp[j]

        col = col + 1
        
    return m


def generateBTMatrixPython(factory, training, verbose=False):
    from pysgpp import DataVector
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    alpha = DataVector(storage.size())
    temp = DataVector(training.getSize())
    
    # create BT matrix
    m = DataVector(training.getSize(), storage.size())
    
    for i in xrange(storage.size()):
        # apply unit vectors
        temp.setAll(0.0)
        alpha.setAll(0.0)
        alpha[i] = 1.0
        b.multTranspose(alpha, training, temp)
        
        #Sets the column in m       
        m.setColumn(i, temp)
        
    return m


def generateBBTMatrixPython(factory, training, verbose=False):
    storage = factory.getStorage()
       
    b = factory.createOperationB()
    
    alpha = DataVector(storage.size())
    erg = DataVector(alpha.getSize())
    temp = DataVector(training.getSize())
    
    # create B matrix
    m = DataVector(storage.size(), storage.size())
    for i in xrange(storage.size()):
        # apply unit vectors
        temp.setAll(0.0)
        erg.setAll(0.0)
        alpha.setAll(0.0)
        alpha[i] = 1.0
        b.multTranspose(alpha, training, temp)
        b.mult(temp, training, erg)
        #Sets the column in m
        m.setColumn(i, erg)
        
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
#    
#    maplematrix = maplematrix + "]):"
#        
#    return maplematrix

def readReferenceMatrix(storage, filename):
    from pysgpp import DataVector
    # read reference matrix
    try:
        fd = gzOpen(filename, 'r')
    except IOError, e:
        fd = None
        
    if not fd:
        fd = gzOpen('../tests/' + filename, 'r')
        
    dat = fd.read().strip()
    fd.close()
    dat = dat.split('\n')
    dat = map(lambda l: l.strip().split(None), dat)

    #print len(dat)
    #print len(dat[0])
    m_ref = DataVector(len(dat), len(dat[0]))
    for i in xrange(len(dat)):
        for j in xrange(len(dat[0])):
            #print float(dat[i][j])
            m_ref[i*len(dat[0]) + j] = float(dat[i][j])

    return m_ref


def readDataVector(filename):
    
    try:
        fin = tools.gzOpen(filename, 'r')
    except IOError, e:
        fin = None
        
    if not fin:
        fin = tools.gzOpen('tests/' + filename, 'r')
    
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

##
# Compares, if two BBT matrices are "almost" equal.
# Has to handle the problem that the underlying grid was ordered
# differently. Uses heuristics, e.g. whether the diagonal elements
# and row and column sums match.
def compareBBTMatrices(m1, m2):
    # check dimensions
 
    n = m1.getSize()

    # check diagonal
    values = []
    for i in range(n):
        values.append(m1[i*n + i])
    values.sort()
    values_ref = []
    for i in range(n):
        values_ref.append(m2[i*n + i])
    values_ref.sort()

    for i in range(n):
        print values_ref[i], values[i]
 #       testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Diagonal %f != %f" % (values[i], values_ref[i]))

    # check row sum
    v = DataVector(n)
    values = []
    for i in range(n):
        m1.getRow(i,v)
        values.append(v.sum())
    values.sort()
    values_ref = []
    for i in range(n):
        m2.getRow(i,v)
        values_ref.append(v.sum())
    values_ref.sort()
    for i in range(n):
        print values_ref[i], values[i]
        #testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Row sum %f != %f" % (values[i], values_ref[i]))

    # check col sum
    v = DataVector(n)
    values = []
    for i in range(n):
        m1.getColumn(i,v)
        values.append(v.sum())
    values.sort()
    values_ref = []
    for i in range(n):
        m2.getColumn(i,v)
        values_ref.append(v.sum())
    values_ref.sort()
    for i in range(n):
        print values_ref[i], values[i]
        #testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Col sum %f != %f" % (values[i], values_ref[i]))


##
# Compares, if two BT matrices are "almost" equal.
# Has to handle the problem that the underlying grid was ordered
# differently. Uses heuristics, e.g. whether the 
# row and column sums match.
def compareBTMatrices(m1, m2):
    # check dimensions
 
    print m1
    print m2
 
    n = m1.getSize() # lines
    m = m1.getDim()  # columns

    # check row sum
    v = DataVector(m)
    values = []
    for i in range(n):
        m1.getRow(i,v)
        values.append(v.sum())
    values.sort()
    values_ref = []
    for i in range(n):
        m2.getRow(i,v)
        values_ref.append(v.sum())
    values_ref.sort()
    for i in range(n):
        print values_ref[i], values[i]
        #testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Row sum %f != %f" % (values[i], values_ref[i]))

    # check col sum
    v = DataVector(n)
    values = []
    for i in range(m):
        m1.getColumn(i,v)
        values.append(v.sum())
    values.sort()
    values_ref = []
    for i in range(m):
        m2.getColumn(i,v)
        values_ref.append(v.sum())
    values_ref.sort()
    for i in range(m):
        print values_ref[i], values[i]
        #testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Col sum %f != %f" % (values[i], values_ref[i]))


def writeMatrixToFile(filename, matrix, n, m):
    fout = file(filename, "w")
        
    for i in range(n):
        for j in range(m):
            fout.write(str(matrix[i*m+j]) + " " )
            
        fout.write("\n")
        
    fout.close()
    
    return

def build_DM_Matrices():
    factory = Grid.createLinearGrid(6)
    level = 3
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    #training = buildTrainingVector(openFile('../datasets/twospirals/twospirals.wieland.arff.gz'))
    #training = buildTrainingVector(openFile('../datasets/ripley/ripleyGarcke.train.arff.gz'))
    training = buildTrainingVector(openFile('../datasets/bupa_liver/liver-disorders_normalized.arff.gz'))
    #training = buildTrainingVector(openFile('../tests/data/data_dim_1_nops_8_float.arff.gz'))
       
    aem = 325
    lam = 0.001
    
    # comparison of matrices
    #m = generateBTMatrixPython(factory, training)
    #m_ref = readReferenceMatrix(factory.getStorage(), 'data/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz')
    #compareBTMatrices(m, m_ref) 
    #writeMatrixToFile('BT_trapezrand_dim_1_17.dat', m, m.getSize(), m.getDim())
    
    print "generating laplacian matrix..."
    C = generateCMatrix(factory, level)
    print C
    print "generating B*B^T matrix..."
    B_res = generateBBTMatrix(factory, training)
    print B_res
    print "multiplying aem*lambda*C..."
    C *= aem
    C *= lam
    print "adding C and B_res..."
    C += B_res
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

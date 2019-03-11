# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest, tools
from pysgpp import createOperationMultipleEvalVectorized

#-------------------------------------------------------------------------------
## Builds the training data vector
# 
# @param data a list of lists that contains the points a the training data set, coordinate-wise
# @return a instance of a DataVector that stores the training data
def buildTrainingVector(data):
    from pysgpp import DataMatrix
    dim = len(data["data"])
    training = DataMatrix(len(data["data"][0]), dim)
    
    # i iterates over the data points, d over the dimension of one data point
    for i in range(len(data["data"][0])):
        for d in range(dim):
            training.set(i, d, data["data"][d][i])
    
    return training


def openFile(filename):
    try:
        data = tools.readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        
    if ("classes" in data) == False:
        print ("No classes found in the given File " + filename + "!")
        
    return data


def generateBBTAVXMatrix(factory, training, verbose=False):
    from pysgpp import DataVector, DataMatrix
    storage = factory.getStorage()
       
    b = createOperationMultipleEvalVectorized(factory,"AVX", training)

    alpha = DataVector(storage.getSize())
    erg = DataVector(len(alpha))
    
    #padding
    remainder = training.getNrows() % 24
    loopCount = 24 - remainder
    if loopCount != 24:
        lastRow = DataVector(training.getNcols())
        for i in range(0, loopCount):
            training.getRow(training.getNrows()-1, lastRow)
            training.resize(training.getNrows()+1)
            training.setRow(training.getNrows()-1, lastRow)

    training.transpose()
    
    temp = DataVector(training.getNcols())
    # create B matrix
    m = DataMatrix(storage.getSize(), storage.getSize())
    for i in range(storage.getSize()):
        # apply unit vectors
        temp.setAll(0.0)
        erg.setAll(0.0)
        alpha.setAll(0.0)
        alpha[i] = 1.0
        b.multVectorized(alpha, temp)
        b.multTransposeVectorized(temp, erg)
        #Sets the column in m
        m.setColumn(i, erg)
        
    return m


def readReferenceMatrix(self, storage, filename):
    from pysgpp import DataVector, DataMatrix
    # read reference matrix
    try:
        fd = tools.gzOpen(filename, 'r')
    except IOError as e:
        fd = None
        
    if not fd:
        fd = tools.gzOpen('tests/' + filename, 'r')
        
    dat = fd.read().strip()
    fd.close()
    dat = dat.split('\n')
    dat = [l.strip().split(None) for l in dat]

    # right number of entries?
    self.assertEqual(storage.getSize(), len(dat))
    self.assertEqual(storage.getSize(), len(dat[0]))

    m_ref = DataMatrix(len(dat), len(dat[0]))
    for i in range(len(dat)):
        for j in range(len(dat[0])):
            m_ref.set(i, j, float(dat[i][j]))

    return m_ref

def readDataVector(filename):
    from pysgpp import DataVector
    
    try:
        fin = tools.gzOpen(filename, 'r')
    except IOError as e:
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
        for i in range(len(values)):
            data[i].append(float(values[i]))
            
    # cleaning up and return
    fin.close()
    return {"data":data, "classes":classes, "filename":filename}

##
# Compares, if two matrices are "almost" equal.
# Has to handle the problem that the underlying grid was ordered
# differently. Uses heuristics, e.g. whether the diagonal elements
# and row and column sums match.
def compareBBTAVXMatrices(testCaseClass, m1, m2):
    from pysgpp import DataVector

    # check dimensions
    testCaseClass.assertEqual(m1.getNrows(), m1.getNcols())
    testCaseClass.assertEqual(m1.getNrows(), m2.getNrows())
    testCaseClass.assertEqual(m1.getSize(), m2.getSize())

    n = m1.getNrows()

    # check diagonal
    values = []
    for i in range(n):
        values.append(m1.get(i,i))
    values.sort()
    values_ref = []
    for i in range(n):
        values_ref.append(m2.get(i,i))
    values_ref.sort()
    for i in range(n):
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], 5, msg="Diagonal %f != %f" % (values[i], values_ref[i]))

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
        #print values_ref[i], values[i]
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], 5, msg="Row sum %f != %f" % (values[i], values_ref[i]))

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
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], 5, msg="Col sum %f != %f" % (values[i], values_ref[i]))

class TestOperationBBTAVXLinear(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(1)
        training = buildTrainingVector(readDataVector('data/data_dim_1_nops_8_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBBTAVXMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/BBT_phi_li_hut_dim_1_nopsgrid_7_float.dat.gz')

        # compare
        compareBBTAVXMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(1)
        training = buildTrainingVector(readDataVector('data/data_dim_1_nops_8_float.arff.gz'))
        level = 5
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBBTAVXMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/BBT_phi_li_hut_dim_1_nopsgrid_31_float.dat.gz')

        # compare
        compareBBTAVXMatrices(self, m, m_ref) 

                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(3)
        training = buildTrainingVector(readDataVector('data/data_dim_3_nops_512_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBBTAVXMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/BBT_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareBBTAVXMatrices(self, m, m_ref) 
  
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(3)
        training = buildTrainingVector(readDataVector('data/data_dim_3_nops_512_float.arff.gz'))
        level = 4
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBBTAVXMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/BBT_phi_li_hut_dim_3_nopsgrid_111_float.dat.gz')

        # compare
        compareBBTAVXMatrices(self, m, m_ref) 
      
                                       
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()

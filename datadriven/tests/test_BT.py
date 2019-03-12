# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import tools
from pysgpp import createOperationMultipleEval

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
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            training.set(i, d, data["data"][d][i])
    
    return training


def openFile(filename):
    try:
        data = tools.readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        
    if data.has_key("classes") == False:
        print ("No classes found in the given File " + filename + "!")
        
    return data


def generateBTMatrix(factory, training, verbose=False):
    from pysgpp import DataVector, DataMatrix
    storage = factory.getStorage()
       
    b = createOperationMultipleEval(factory, training)
    
    alpha = DataVector(storage.getSize())
    temp = DataVector(training.getNrows())
    
    # create BT matrix
    m = DataMatrix(training.getNrows(), storage.getSize())
    
    for i in xrange(storage.getSize()):
      
        # apply unit vectors
        temp.setAll(0.0)
        alpha.setAll(0.0)
        alpha[i] = 1.0
        b.mult(alpha, temp)
        
        #Sets the column in m       
        m.setColumn(i, temp)
        
    return m


def readReferenceMatrix(self, storage, filename):
    from pysgpp import DataMatrix
    # read reference matrix
    try:
        fd = tools.gzOpen(filename, 'r')
    except IOError, e:
        fd = None
        
    if not fd:
        fd = tools.gzOpen('tests/' + filename, 'r')
        
    dat = fd.read().strip()
    fd.close()
    dat = dat.split('\n')
    dat = map(lambda l: l.strip().split(None), dat)

    # right number of entries?
    self.assertEqual(storage.getSize(), len(dat[0]))

    m_ref = DataMatrix(len(dat), len(dat[0]))
    for i in xrange(len(dat)):
        for j in xrange(len(dat[0])):
            m_ref.set(i, j, float(dat[i][j]))

    return m_ref

def readDataVector(filename):
    from pysgpp import DataVector
    
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
# Compares, if two matrices are "almost" equal.
# Has to handle the problem that the underlying grid was ordered
# differently. Uses heuristics, e.g. whether the diagonal elements
# and row and column sums match.
def compareBTMatrices(testCaseClass, m1, m2):
    from pysgpp import DataVector

    places = 5

    # check dimensions
    testCaseClass.assertEqual(m1.getNrows(), m2.getNrows())
    testCaseClass.assertEqual(m1.getNcols(), m2.getNcols())

    n = m1.getNrows() # lines
    m = m1.getNcols() # columns

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
        #print values_ref[i], values[i]
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], places=places, msg="Row sum %f != %f" % (values[i], values_ref[i]))

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
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], places=places, msg="Col sum %f != %f" % (values[i], values_ref[i]))



class TestOperationBTModLinear(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(1)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_ausgeklappt_dim_1_nopsgrid_7_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(1)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 5
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_ausgeklappt_dim_1_nopsgrid_31_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(3)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_ausgeklappt_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 
  
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(3)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 4
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_ausgeklappt_dim_3_nopsgrid_111_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 
             

class TestOperationBTLinear(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(1)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_dim_1_nopsgrid_7_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(1)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 5
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_dim_1_nopsgrid_31_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(3)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 
  
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(3)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 4
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_dim_3_nopsgrid_111_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 
             

class TestOperationBTLinearBoundary(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1, 0)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 4
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1, 0)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 5
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3, 0)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 
  
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3, 0)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 4
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref)     
        

class TestOperationBTLinearTruncatedBoundary(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 4
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_1_nops_8_float.arff.gz'))
        level = 5
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 

                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 2
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref) 
  
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3)
        training = buildTrainingVector(readDataVector('../datasets/BT_BBT/data_dim_3_nops_512_float.arff.gz'))
        level = 3
        gen = factory.getGenerator()
        gen.regular(level)

        m = generateBTMatrix(factory, training)
        m_ref = readReferenceMatrix(self, factory.getStorage(), '../datasets/BT_BBT/BT_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz')

        # compare
        compareBTMatrices(self, m, m_ref)  
        
                                       
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()

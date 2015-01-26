// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)####################################################################

import unittest, tools
from pysgpp import *

def generateLaplaceMatrix(factory, level, verbose=False):
    from pysgpp import DataVector, DataMatrix
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    laplace = createOperationLaplace(factory)
    
    # create vector
    alpha = DataVector(storage.size())
    erg = DataVector(storage.size())

    # create stiffness matrix
    m = DataMatrix(storage.size(), storage.size())
    m.setAll(0)
    for i in xrange(storage.size()):
        # apply unit vectors
        alpha.setAll(0)
        alpha[i] = 1.0
        laplace.mult(alpha, erg)
        if verbose:
            print erg, erg.sum()
        m.setColumn(i, erg)

    return m

def generateLaplaceEnhancedMatrix(factory, level, verbose=False):
    from pysgpp import DataVector, DataMatrix
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(level)
    
    laplace = createOperationLaplaceEnhanced(factory)
    
    # create vector
    alpha = DataVector(storage.size())
    erg = DataVector(storage.size())

    # create stiffness matrix
    m = DataMatrix(storage.size(), storage.size())
    m.setAll(0)
    for i in xrange(storage.size()):
        # apply unit vectors
        alpha.setAll(0)
        alpha[i] = 1.0
        laplace.mult(alpha, erg)
        if verbose:
            print erg, erg.sum()
        m.setColumn(i, erg)

    return m

# is currently only available for MPI
#####################################
#def generateLaplaceVectorizedMatrix(factory, level, verbose=False):
#    from pysgpp import DataVector, DataMatrix
#    storage = factory.getStorage()
#    
#    gen = factory.createGridGenerator()
#    gen.regular(level)
#    
#    laplace = createOperationLaplaceVectorized(factory)
#    
#    # create vector
#    alpha = DataVector(storage.size())
#    erg = DataVector(storage.size())
#
#    # create stiffness matrix
#    m = DataMatrix(storage.size(), storage.size())
#    m.setAll(0)
#    for i in xrange(storage.size()):
#        # apply unit vectors
#        alpha.setAll(0)
#        alpha[i] = 1.0
#        laplace.mult(alpha, erg)
#        if verbose:
#            print erg, erg.sum()
#        m.setColumn(i, erg)
#
#    return m

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
    self.assertEqual(storage.size(), len(dat))
    self.assertEqual(storage.size(), len(dat[0]))

    m_ref = DataMatrix(len(dat), len(dat[0]))
    for i in xrange(len(dat)):
        for j in xrange(len(dat[0])):
            m_ref.set(i, j, float(dat[i][j]))
    
    return m_ref

##
# Compares, if two stiffness matrices are "almost" equal.
# Has to handle the problem that the underlying grid was ordered
# differently. Uses heuristics, e.g. whether the diagonal elements
# and row and column sums match.
def compareStiffnessMatrices(testCaseClass, m1, m2):
    from pysgpp import DataVector, DataMatrix

    # check dimensions
    testCaseClass.assertEqual(m1.getNrows(), m1.getNcols())
    testCaseClass.assertEqual(m1.getNrows(), m2.getNrows())
    testCaseClass.assertEqual(m1.getNcols(), m2.getNcols())

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
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Diagonal %f != %f" % (values[i], values_ref[i]))

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
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Row sum %f != %f" % (values[i], values_ref[i]))

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
        testCaseClass.assertAlmostEqual(values[i], values_ref[i], msg="Col sum %f != %f in Col %f" % (values[i], values_ref[i], i))



class TestOperationLaplaceLinear(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D(self):
        from pysgpp import Grid, DataVector, DataMatrix
        
        factory = Grid.createLinearGrid(1)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(7)
        
        laplace = createOperationLaplace(factory)

        
        alpha = DataVector(storage.size())
        result = DataVector(storage.size())
        
        alpha.setAll(1.0)
        
        laplace.mult(alpha, result)
        
        for seq in xrange(storage.size()):
            index = storage.get(seq)
            level, _ = index.get(0)
            self.failUnlessAlmostEqual(result[seq], pow(2.0, level+1))

        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD(self):
        
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(3)

        m = generateLaplaceMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)

class TestOperationLaplaceEnhancedLinear(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D(self):
        from pysgpp import Grid, DataVector, DataMatrix
        
        factory = Grid.createLinearGrid(1)
        storage = factory.getStorage()
        
        gen = factory.createGridGenerator()
        gen.regular(7)
        
        laplace = createOperationLaplaceEnhanced(factory)
        
        alpha = DataVector(storage.size())
        result = DataVector(storage.size())
        
        alpha.setAll(1.0)
        
        laplace.mult(alpha, result)
        
        for seq in xrange(storage.size()):
            index = storage.get(seq)
            level, _ = index.get(0)
            self.failUnlessAlmostEqual(result[seq], pow(2.0, level+1))

        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD(self):
        
        from pysgpp import Grid
        
        factory = Grid.createLinearGrid(3)

        m = generateLaplaceEnhancedMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
        
#class TestOperationLaplaceLinearVectorized(unittest.TestCase):
#    ##
#    # Test laplace for regular sparse grid in 1d using linear hat functions
#    def testHatRegular1D(self):
#        from pysgpp import Grid, DataVector, DataMatrix
#        
#        factory = Grid.createLinearGrid(1)
#        storage = factory.getStorage()
#        
#        gen = factory.createGridGenerator()
#        gen.regular(7)
#        
#        laplace = createOperationLaplaceVectorized(factory)
#      
#        
#        alpha = DataVector(storage.size())
#        result = DataVector(storage.size())
#        
#        alpha.setAll(1.0)
#        
#        laplace.mult(alpha, result)
#        
#        for seq in xrange(storage.size()):
#            index = storage.get(seq)
#            level, _ = index.get(0)
#            self.failUnlessAlmostEqual(result[seq], pow(2.0, level+1))
#
#        
#    ##
#    # Test regular sparse grid dD, normal hat basis functions.
#    def testHatRegulardD(self):
#        
#        from pysgpp import Grid
#        
#        factory = Grid.createLinearGrid(3)
#
#        m = generateLaplaceVectorizedMatrix(factory, 3)
#        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz')
#
#        # compare
#        compareStiffnessMatrices(self, m, m_ref)
#
   
class TestOperationLaplaceModLinear(unittest.TestCase):
    def testHatRegular1D(self):
        from pysgpp import Grid, DataVector, DataMatrix
        
        factory = Grid.createModLinearGrid(1)

        m = generateLaplaceMatrix(factory, 5)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_ausgeklappt_dim_1_nopsgrid_31_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
   
    def testHatRegulardD(self):
        from pysgpp import Grid
        
        factory = Grid.createModLinearGrid(3)
        #print "------------------------------------------------"
        m = generateLaplaceMatrix(factory, 3)
        #print "------------------------------------------------"
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_ausgeklappt_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
        
        
class TestOperationLaplaceLinearTrapezoidBoundary(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(1)

        m = generateLaplaceMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(1)

        m = generateLaplaceMatrix(factory, 5)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 
        
                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(3)

        m = generateLaplaceMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)  
           
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(3)

        m = generateLaplaceMatrix(factory, 2)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
        
    
class TestOperationLaplaceEnhancedLinearTrapezoidBoundary(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(1)

        m = generateLaplaceEnhancedMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(1)

        m = generateLaplaceEnhancedMatrix(factory, 5)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 
        
                
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):  
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(3)

        m = generateLaplaceEnhancedMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)  
           
        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearTrapezoidBoundaryGrid(3)

        m = generateLaplaceEnhancedMatrix(factory, 2)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)

                
class TestOperationLaplaceLinearBoundary(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1)

        m = generateLaplaceMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1)

        m = generateLaplaceMatrix(factory, 5)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 

        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3)

        m = generateLaplaceMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)  
        
    
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3)

        m = generateLaplaceMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)    
        
        
class TestOperationLaplaceEnhancedLinearBoundary(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1)

        m = generateLaplaceEnhancedMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 

  
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testHatRegular1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(1)

        m = generateLaplaceEnhancedMatrix(factory, 5)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 

        
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_one(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3)

        m = generateLaplaceEnhancedMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)  
        
    
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHatRegulardD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createLinearBoundaryGrid(3)

        m = generateLaplaceEnhancedMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)    


class TestOperationLaplacePrewavelet(unittest.TestCase):
    ##
    # Test laplace for regular sparse grid in 1d using linear hat functions
    def testPrewavelet1D_one(self):
        from pysgpp import Grid
        
        factory = Grid.createPrewaveletGrid(1)

        m = generateLaplaceMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_prewavelet_dim_1_nopsgrid_15_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref) 
        
    def testPrewavelet1D_two(self):
        from pysgpp import Grid
        
        factory = Grid.createPrewaveletGrid(1)

        m = generateLaplaceMatrix(factory, 5)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_prewavelet_dim_1_nopsgrid_31_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
        
    def testPrewaveletdD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createPrewaveletGrid(3)

        m = generateLaplaceMatrix(factory, 3)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_prewavelet_dim_3_nopsgrid_31_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
        
    def testPrewaveletdD_two(self):
        from pysgpp import Grid
        
        factory = Grid.createPrewaveletGrid(3)

        m = generateLaplaceMatrix(factory, 4)
        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_prewavelet_dim_3_nopsgrid_111_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
        
    def testPrewaveletAdaptivedD_two(self):
        from pysgpp import Grid, DataVector, SurplusRefinementFunctor, DataMatrix
        
        factory = Grid.createPrewaveletGrid(4)
        level = 2
        gen = factory.createGridGenerator()
        gen.regular(level)
        
        alpha = DataVector(factory.getStorage().size())
        for i in xrange(factory.getStorage().size()):
            alpha[i]=i+1
        gen.refine(SurplusRefinementFunctor(alpha,1));


        storage = factory.getStorage()
    
        laplace = createOperationLaplace(factory)
    
    # create vector
        alpha = DataVector(storage.size())
        erg = DataVector(storage.size())

    # create stiffness matrix
        m = DataMatrix(storage.size(), storage.size())
        m.setAll(0)
        for i in xrange(storage.size()):
        # apply unit vectors
            alpha.setAll(0)
            alpha[i] = 1.0
            laplace.mult(alpha, erg)
            m.setColumn(i, erg)

        m_ref = readReferenceMatrix(self, factory.getStorage(), 'data/C_laplace_prewavelet_dim_4_nopsgrid_17_adapt_float.dat.gz')

        # compare
        compareStiffnessMatrices(self, m, m_ref)
                 
        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main(verbosity=9)


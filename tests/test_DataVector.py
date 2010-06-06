###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de)


import unittest

##
# @package tests.test_DataVector
# Contains class test_DataVector::TestDataVector with unittests for DataVector

##
# Class with unittests for DataVector
#
# @ingroup tests
#
# @test Unittests for DataVector
# @todo (pflueged) Test cases for remaining functionality of DataVector
class TestDataVector(unittest.TestCase):

    ## 
    # Set up, create random DataVector and corresponding Python data structures.
    # @test DataVector::get(), DataVector::set()
    def setUp(self):
        from pysgpp import DataVector
        import random

        ## number of rows
        self.nrows = 5
        ## number of columns
        self.ncols = 4
        ## number of entries
        self.N = self.nrows*self.ncols
        ## random list of lists
        self.l_rand = [[2*(random.random()-0.5) for j in xrange(self.ncols)] for i in xrange(self.nrows)]
        ## same as l_rand, but flattened
        self.l_rand_total = []
        for li in self.l_rand:
            self.l_rand_total.extend(li)
#        ## Data Vector, corresponding to l_rand
#        self.d_rand = DataVector(self.nrows,self.ncols)
#        for i in xrange(self.N):
#            self.d_rand[i] = self.l_rand_total[i]
#
#        for i in xrange(self.N):
#            self.assertEqual(self.d_rand[i], self.l_rand_total[i])
        ## Data Vector, corresponding to l_rand
        self.d_rand = DataVector(self.N)
        for i in xrange(self.N):
            self.d_rand[i] = self.l_rand_total[i]

        for i in xrange(self.N):
            self.assertEqual(self.d_rand[i], self.l_rand_total[i])

    ##
    # Constructors4.
    # @test DataVector::DataVector(size_t size), DataVector::DataVector(size_t size, size_t dim), DataVector::DataVector(DataVectorDefinition &DataVectorDef), DataVector::getSize(), DataVector::getDim(), DataVector::getSize()
    # @todo (pflueged) DataVector::DataVector(double *input, size_t size, size_t dim)
    def testConstructor(self):
        from pysgpp import DataVector
        
        d = DataVector(2)
        self.assertEqual(len(d), 2) # getSize()
        
#        d = DataVector(2,3)
#        self.assertEqual(d.getSize(), 2)
#        self.assertEqual(d.getDim(), 3)
#        self.assertEqual(len(d), 2*3) # getSize()
#
#        d2 = DataVector(self.d_rand)
#        for i in xrange(self.N):
#            self.assertEqual(d2[i], self.d_rand[i])
#        self.assertEqual(d2.getSize(), self.nrows)
#        self.assertEqual(d2.getDim(), self.ncols)
#        self.assertEqual(len(d2), self.N)
#        d2[self.ncols] = -4.0
#        self.assertNotEqual(d2[self.ncols], self.d_rand[self.ncols])

    ##
    # Min, Max operations.
    # @test DataVector::min(int d), DataVector::max(int d), DataVector::minmax(int d, double *min, double *max), DataVector::min(), DataVector::max()
    def testMinMax(self):

#        # test dimension-dependent min, max
#        for j in xrange(self.ncols):
#            minj = min([self.l_rand[i][j] for i in xrange(self.nrows)])
#            maxj = max([self.l_rand[i][j] for i in xrange(self.nrows)])
#            self.assertEqual(self.d_rand.min(j), minj)
#            self.assertEqual(self.d_rand.max(j), maxj)
#            mi, ma = self.d_rand.minmax(j)
#            self.assertEqual(mi, minj)
#            self.assertEqual(ma, maxj)

        # test global min, max
        self.assertEqual(self.d_rand.min(), min(self.l_rand_total))
        self.assertEqual(self.d_rand.max(), max(self.l_rand_total))
   

    ##
    # Operations on DataVectors.
    # @test DataVector::sum(), DataVector::sqr(), DataVector::abs(), DataVector::componentwise_mult(), DataVector::componentwise_div()
    def testOps(self):
        from pysgpp import DataVector
        # sum
        self.assertAlmostEqual(self.d_rand.sum(), sum(self.l_rand_total))

        # sqr
        d = DataVector(self.d_rand)
        d.sqr()
        for i in xrange(self.N):
            self.assertEqual(self.d_rand[i]**2, d[i])

        # abs
        d = DataVector(self.d_rand)
        d.abs()
        for i in xrange(self.N):
            self.assertEqual(abs(self.d_rand[i]), d[i])

        # componentwise_mult
        d = DataVector(self.d_rand)
#	d2 = DataVector(self.nrows, self.ncols)
	d2 = DataVector(self.N)
        for i in xrange(self.N):
            d2[i] = i
	d.componentwise_mult(d2)
        for i in xrange(self.N):
            self.assertEqual(self.d_rand[i]*i, d[i])

        # componentwise_div
        d = DataVector(self.d_rand)
        for i in xrange(self.N):
            d2[i] = i+1
	d.componentwise_div(d2)
        for i in xrange(self.N):
            self.assertEqual(self.d_rand[i]/(i+1), d[i])

    ##
    # Vector-Operations
    # @test DataVector::dotProduct(DataVector &vec)
    def testDotProduct(self):
        from pysgpp import DataVector
        
        x = 0
        
        d = DataVector(3)
        for i in xrange(len(d)):
            d[i] = i + 1
            x += d[i] * d[i]
            
        self.assertEqual(d.dotProduct(d), x)

# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()

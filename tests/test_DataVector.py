#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007-2009 Dirk Pflueger (pflueged@in.tum.de)                #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################


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
    # @test DataVector::get()
    # @test DataVector::set()
    def setUp(self):
        from pysgpp import DataVector
        import random

        self.nrows = 5
        self.ncols = 4
        self.N = self.nrows*self.ncols
        # random list
        self.l_rand = [[2*(random.random()-0.5) for j in xrange(self.ncols)] for i in xrange(self.nrows)]
        self.l_rand_total = []
        for li in self.l_rand:
            self.l_rand_total.extend(li)
        # Data Vector
        self.d_rand = DataVector(self.nrows,self.ncols)
        for i in xrange(self.N):
            self.d_rand[i] = self.l_rand_total[i]

        for i in xrange(self.N):
            self.assertEqual(self.d_rand[i], self.l_rand_total[i])

    ##
    # Constructors.
    # @test DataVector::DataVector(size_t size)
    # @test DataVector::DataVector(size_t size, size_t dim)
    # @todo (pflueged) DataVector::DataVector(double *input, size_t size, size_t dim), DataVector::DataVector(DataVectorDefinition &DataVectorDef)
    def testConstructor(self):
        from pysgpp import DataVector
        
        d = DataVector(2)
        self.assertEqual(d.getSize(), 2)
        
        d = DataVector(2,2)
        self.assertEqual(d.getSize(), 2)
        self.assertEqual(d.getDim(), 2)
        self.assertEqual(len(d), 4)

    ##
    # Min, Max operations.
    # @test DataVector::min(int d)
    # @test DataVector::max(int d)
    # @test DataVector::minmax(int d, double *min, double *max)
    # @test DataVector::min()
    # @test DataVector::max()
    def testMinMax(self):

        # test dimension-dependent min, max
        for j in xrange(self.ncols):
            minj = min([self.l_rand[i][j] for i in xrange(self.nrows)])
            maxj = max([self.l_rand[i][j] for i in xrange(self.nrows)])
            self.assertEqual(self.d_rand.min(j), minj)
            self.assertEqual(self.d_rand.max(j), maxj)
            mi, ma = self.d_rand.minmax(j)
            self.assertEqual(mi, minj)
            self.assertEqual(ma, maxj)

        # test global min, max
        self.assertEqual(self.d_rand.min(), min(self.l_rand_total))
        self.assertEqual(self.d_rand.max(), max(self.l_rand_total))
   

    ##
    # Operations on DataVectors.
    # @test DataVector::sum()
    # @test DataVector::sqr()
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

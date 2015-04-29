# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp
import random
import tempfile

class TestTools(unittest.TestCase):
    def setUp(self):
        """Initialize the test case."""
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def gridEqualityTest(self, grid1, grid2):
        grid1Storage = grid1.getStorage()
        grid2Storage = grid2.getStorage()
        d = grid1Storage.dim()
        self.assertEqual(grid1Storage.size(), grid2Storage.size())
        self.assertEqual(d, grid2Storage.dim())
        
        for k in range(grid1Storage.size()):
            for t in range(d):
                self.assertEqual(grid1Storage.get(k).getLevel(t),
                                 grid2Storage.get(k).getLevel(t))
                self.assertEqual(grid1Storage.get(k).getIndex(t),
                                 grid2Storage.get(k).getIndex(t))
    
    def testFileIOReadWriteGrid(self):
        """Test SGPP::optimization::file_io::readGrid/writeGrid."""
        random.seed(42)
        
        for d in range(1, 6):
            grid1 = pysgpp.Grid.createLinearGrid(d)
            gridGen = grid1.createGridGenerator()
            gridGen.regular(3)
            
            with tempfile.NamedTemporaryFile() as f:
                pysgpp.OptFileIOWriteGrid(f.name, grid1.getStorage())
                grid2 = pysgpp.Grid.createLinearGrid(d)
                pysgpp.OptFileIOReadGrid(f.name, grid2.getStorage())
            
            self.gridEqualityTest(grid1, grid2)
            
            functionValues1 = pysgpp.DataVector(
                    [random.uniform(0.0, 1.0)
                     for k in range(grid1.getStorage().size())])
            
            with tempfile.NamedTemporaryFile() as f:
                pysgpp.OptFileIOWriteGrid(f.name, grid1.getStorage(),
                                          functionValues1)
                grid2 = pysgpp.Grid.createLinearGrid(d)
                functionValues2 = pysgpp.DataVector(0)
                pysgpp.OptFileIOReadGrid(f.name, grid2.getStorage(),
                                         functionValues2)
            
            self.gridEqualityTest(grid1, grid2)
            self.assertEqual(len(functionValues1), len(functionValues2))
            for k in range(len(functionValues1)):
                if pysgpp.cvar.USING_DOUBLE_PRECISION:
                    self.assertEqual(functionValues1[k], functionValues2[k])
                else:
                    self.assertAlmostEqual(functionValues1[k], functionValues2[k])
    
    def testFileIOReadWriteMatrix(self):
        """Test SGPP::optimization::file_io::readMatrix/writeMatrix."""
        random.seed(42)
        m1, n1 = (100, 200)
        
        # test read/write with std::vector<SGPP::float_t>
        A1 = pysgpp.DoubleVector()
        for i in range(m1):
            for j in range(n1):
                A1.push_back(random.uniform(0.0, 1.0))
        
        with tempfile.NamedTemporaryFile() as f:
            pysgpp.OptFileIOWriteMatrix(f.name, A1, m1, n1)
            A2 = pysgpp.DoubleVector()
            m2, n2 = pysgpp.OptFileIOReadMatrix(f.name, A2)
        
        self.assertEqual(m1, m2)
        self.assertEqual(n1, n2)
        self.assertEqual(len(A1), len(A2))
        for k in range(len(A1)):
            self.assertEqual(A1[k], A2[k])
        
        # test read/write with DataMatrix
        A1 = pysgpp.DataMatrix(m1, n1)
        for i in range(m1):
            for j in range(n1):
                A1.set(i, j, random.uniform(0.0, 1.0))
        
        with tempfile.NamedTemporaryFile() as f:
            pysgpp.OptFileIOWriteMatrix(f.name, A1)
            A2 = pysgpp.DataMatrix(0, 0)
            pysgpp.OptFileIOReadMatrix(f.name, A2)
        
        self.assertEqual(m1, m2)
        self.assertEqual(n1, n2)
        self.assertEqual(A1.getNrows(), A2.getNrows())
        self.assertEqual(A1.getNcols(), A2.getNcols())
        for i in range(m1):
            for j in range(n1):
                self.assertEqual(A1.get(i, j), A2.get(i, j))
    
    def testFileIOReadWriteVector(self):
        """Test SGPP::optimization::file_io::readVector/writeVector."""
        random.seed(42)
        n = 100
        
        # test read/write with std::vector<SGPP::float_t>
        v1 = pysgpp.DoubleVector()
        for i in range(n):
            v1.push_back(random.uniform(0.0, 1.0))
        
        with tempfile.NamedTemporaryFile() as f:
            pysgpp.OptFileIOWriteVector(f.name, v1)
            v2 = pysgpp.DoubleVector()
            pysgpp.OptFileIOReadVector(f.name, v2)
        
        self.assertEqual(len(v1), len(v2))
        for k in range(len(v1)):
            self.assertEqual(v1[k], v2[k])
        
        # test read/write with DataVector
        v1 = pysgpp.DataVector(n)
        for i in range(n):
            v1[i] = random.uniform(0.0, 1.0)
        
        with tempfile.NamedTemporaryFile() as f:
            pysgpp.OptFileIOWriteVector(f.name, v1)
            v2 = pysgpp.DataVector(0)
            pysgpp.OptFileIOReadVector(f.name, v2)
        
        self.assertEqual(len(v1), len(v2))
        for i in range(n):
            self.assertEqual(v1[i], v2[i])
    
    def testRandomNumberGenerator(self):
        """Test SGPP::optimization::RandomNumberGenerator."""
        seed = 42
        N = 10000
        
        # set and test seed getting/setting
        pysgpp.cvar.OptRNGInstance.setSeed()
        pysgpp.cvar.OptRNGInstance.setSeed(seed)
        self.assertEqual(pysgpp.cvar.OptRNGInstance.getSeed(), seed)
        
        # test continuous uniform random numbers
        numbers = [pysgpp.cvar.OptRNGInstance.getUniformRN()
                   for i in range(N)]
        mean = sum(numbers) / len(numbers)
        var = sum([(number - mean)**2 for number in numbers]) / len(numbers)
        
        self.assertTrue(all([isinstance(number, float) for number in numbers]))
        self.assertTrue(all([0.0 <= number <= 1.0 for number in numbers]))
        self.assertAlmostEqual(mean, 0.5, places=2)
        self.assertAlmostEqual(var, 1.0/12.0, places=2)
        
        # test Gaussian random numbers
        mus = [0.0, 12.3, -42.0, 13.37]
        sigmas = [1.0, 2.6, 8.1, 0.3]
        
        for mu, sigma in zip(mus, sigmas):
            numbers = [pysgpp.cvar.OptRNGInstance.getGaussianRN(sigma, mu)
                       for i in range(N)]
            mean = sum(numbers) / len(numbers)
            var = sum([(number - mean)**2 for number in numbers]) / \
                  len(numbers)
                  
            self.assertTrue(all([isinstance(number, float)
                                 for number in numbers]))
            self.assertAlmostEqual(mean, mu, delta=sigma*0.1)
            self.assertAlmostEqual(var, sigma*sigma, delta=sigma*sigma*0.1)
        
        # test discrete uniform random numbers
        for k in range(1, 11):
            numbers = [pysgpp.cvar.OptRNGInstance.getUniformIndexRN(k)
                       for i in range(N)]
            mean = float(sum(numbers)) / float(len(numbers))
            var = float(sum([(number - mean)**2 for number in numbers])) / \
                  float(len(numbers))
            
            self.assertTrue(all([isinstance(number, int) or
                                 isinstance(number, long)
                                 for number in numbers]))
            self.assertTrue(all([0 <= number <= k-1 for number in numbers]))
            self.assertAlmostEqual(mean, (k-1.0)/2.0, delta=k*0.01)
            self.assertAlmostEqual(var, (k*k-1.0)/12.0, delta=k*k*0.01)
    
    def orthogonalityTest(self, A):
        n = A.getNrows()
        self.assertEqual(n, A.getNcols())
        
        for i in range(n):
            for j in range(i+1):
                entry = 0.0
                for l in range(n):
                    entry += A.get(i, l) * A.get(j, l)
                self.assertAlmostEqual(entry, 1.0 if (i == j) else 0.0)
    
    def symmetryTest(self, A):
        n = A.getNrows()
        self.assertEqual(n, A.getNcols())
        
        for i in range(n):
            for j in range(i):
                self.assertAlmostEqual(A.get(i, j), A.get(j, i))
    
    def similiarityTest(self, A, V, B):
        n = A.getNrows()
        self.assertEqual(n, A.getNcols())
        self.assertEqual(n, V.getNrows())
        self.assertEqual(n, V.getNcols())
        self.assertEqual(n, B.getNrows())
        self.assertEqual(n, B.getNcols())
        
        for i in range(n):
            for j in range(n):
                entry1 = 0.0
                entry2 = 0.0
                for l in range(n):
                    entry1 += A.get(i, l) * V.get(l, j)
                    entry2 += V.get(i, l) * B.get(l, j)
                self.assertAlmostEqual(entry1, entry2)
    
    def generateRandomMatrix(self, n):
        A = pysgpp.DataMatrix(n, n)
        for i in range(n):
            for j in range(n):
                A.set(i, j, random.gauss(0.0, 1.0))
        return A
    
    def testMathHouseholderTransformation(self):
        """Test SGPP::optimization::math::householderTransformation."""
        random.seed(42)
        n = 20
        A = self.generateRandomMatrix(n)
        
        for p, q in [(0, 0), (2, 1), (5, 9)]:
            m = n - p
            Q = pysgpp.DataMatrix(m, m)
            pysgpp.OptMathHouseholderTransformation(A, p, q, Q)
            self.symmetryTest(Q)
            self.orthogonalityTest(Q)
            
            for i in range(p+1, n):
                entry = 0.0
                for l in range(p, n):
                    entry += Q.get(i-p, l-p) * A.get(l, q)
                self.assertAlmostEqual(entry, 0.0)
    
    def testMathHessenbergForm(self):
        """Test SGPP::optimization::math::hessenbergForm."""
        random.seed(42)
        n = 20
        A = self.generateRandomMatrix(n)
        H = pysgpp.DataMatrix(A)
        V = pysgpp.DataMatrix(n, n)
        pysgpp.OptMathHessenbergForm(H, V)
        self.orthogonalityTest(V)
        self.similiarityTest(A, V, H)
        
        for i in range(n):
            for j in range(i-2):
                self.assertAlmostEqual(H.get(i, j), 0.0)
    
    def testMathQRDecomposition(self):
        """Test SGPP::optimization::math::QRDecomposition."""
        random.seed(42)
        n = 20
        A = self.generateRandomMatrix(n)
        R = pysgpp.DataMatrix(A)
        Q = pysgpp.DataMatrix(n, n)
        pysgpp.OptMathQRDecomposition(R, Q)
        self.orthogonalityTest(Q)
        
        for i in range(n):
            for j in range(i-1):
                self.assertAlmostEqual(R.get(i, j), 0.0)
        
        for i in range(n):
            for j in range(n):
                entry = 0.0
                for l in range(n):
                    entry += Q.get(i, l) * R.get(l, j)
                self.assertAlmostEqual(A.get(i, j), entry)
    
    def testMathSchurDecomposition(self):
        """Test SGPP::optimization::math::schurDecomposition."""
        random.seed(42)
        n = 20
        A = self.generateRandomMatrix(n)
        
        for i in range(n):
            for j in range(i):
                A.set(i, j, A.get(j, i))
        
        S = pysgpp.DataMatrix(A)
        V = pysgpp.DataMatrix(n, n)
        pysgpp.OptMathSchurDecomposition(S, V)
        self.orthogonalityTest(V)
        self.similiarityTest(A, V, S)
        
        for i in range(n):
            for j in range(n):
                if i != j:
                    self.assertAlmostEqual(S.get(i, j), 0.0)

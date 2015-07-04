# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp

class TestGridGen(unittest.TestCase):
    def setUp(self):
        """Initialize the test case."""
        # disable status output
        pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
        # disable multi-threading
        pysgpp.omp_set_num_threads(1)
    
    def testIterativeGridGenerators(self):
        """Test SGPP::optimization iterative grid generators."""
        pysgpp.cvar.OptRNGInstance.setSeed(42)
        d = 2
        N = 200
        
        f = pysgpp.OptRosenbrock(d)
        f.generateDisplacement()
        
        # Test All The Grids!
        grids = [pysgpp.Grid.createBsplineGrid(d, 3),
                 pysgpp.Grid.createBsplineTruncatedBoundaryGrid(d, 3),
                 pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
                 pysgpp.Grid.createModBsplineGrid(d, 3),
                 pysgpp.Grid.createLinearGrid(d),
                 pysgpp.Grid.createLinearTruncatedBoundaryGrid(d),
                 pysgpp.Grid.createLinearClenshawCurtisGrid(d),
                 pysgpp.Grid.createModLinearGrid(d),
                 pysgpp.Grid.createWaveletGrid(d),
                 pysgpp.Grid.createWaveletTruncatedBoundaryGrid(d),
                 pysgpp.Grid.createModWaveletGrid(d),
                 pysgpp.Grid.createFundamentalSplineGrid(d, 3),
                 pysgpp.Grid.createModFundamentalSplineGrid(d, 3)]
        
        for grid in grids:
            # repeat for Ritter-Novak and linear surplus grid generation
            for k in range(3):
                # empty grid
                grid.getStorage().emptyStorage()
                if k == 0:
                    gridgen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, 0.85)
                elif k == 1:
                    gridgen = pysgpp.OptIterativeGridGeneratorLinearSurplus(f, grid, N, 0.2)
                else:
                    gridgen = pysgpp.OptIterativeGridGeneratorSOO(f, grid, N, 0.5)
                
                # generate grid
                self.assertTrue(gridgen.generate())
                
                # test grid size
                n = grid.getSize()
                self.assertLessEqual(n, N)
                
                # test size of function value vector
                function_values = gridgen.getFunctionValues()
                self.assertEqual(n, len(function_values))
                
                for i in range(n):
                    x = pysgpp.DataVector([grid.getStorage().get(i).getCoord(t) for t in range(d)])
                    # test function value
                    self.assertAlmostEqual(function_values[i], f.eval(x)) 

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import pysgpp

'''
In this example we create a B-spline based surrogate for a function with 
vectorial output. We do this through the convenience class 
'SplineResponseSurface'. For functions with vectorial output the analogous 
class 'SplineResponseSurfaceVector' can be used.
The advantage of these classes is that only few lines of code suffice to create
the surrogate and to calculate various quantities of interest.
'''

# First we set up an objective function. This must be wrapped for SG++


class objFuncSGpp(pysgpp.ScalarFunction):

    def __init__(self, dim):
        super(objFuncSGpp, self).__init__(dim)

    def eval(self, x):
        # input x: pysgpp.DataVector
        return x[0]*x[1]+x[1]


# create an instance of the objective
dim = 2
objFunc = objFuncSGpp(dim)
# set the objectives domain
lb = pysgpp.DataVector([0, 0])  # domain has lower bounds 0
ub = pysgpp.DataVector([1, 1])  # domain has upper bounds 1

# set up response surface
degree = 3
gridType = 'nakBsplineBoundary'

# create a response surface object
reSurf = pysgpp.SplineResponseSurface(
    objFunc, lb, ub, pysgpp.Grid.stringToGridType(gridType), degree)

# create surrogate with regular sparse grid
# reSurf.regular(2)

# create surrogate with spatially adaptive sparse grid
numPoints = 30      # max number of grid points
initialLevel = 1    # initial level
numRefine = 5       # number of grid points refined in each step
verbose = False     # verbosity of subroutines
reSurf.surplusAdaptive(numPoints, initialLevel, numRefine, verbose)

# Now the surrogate has been created and we can use it to calculate various
# quantities of interest

# evaluate the surrogate
evalPoint = pysgpp.DataVector([0.3, 0.6])
print(f'reSurf: {reSurf.eval(evalPoint)}     truth: {objFunc.eval(evalPoint)}')

# evaluate the surrogate's gradient
gradient = pysgpp.DataVector(dim, 0)
reSurf.evalGradient(evalPoint, gradient)
print(f'gradient: {gradient.toString()}')

# integrate the surrogate
print(f'integral: {reSurf.getIntegral()}')

# calculate stochastic moments
# first specify distibutions
pdfs = pysgpp.DistributionsVector()
pdfs.push_back(pysgpp.DistributionNormal(0.5, 0.1))
pdfs.push_back(pysgpp.DistributionUniform(0, 1))
# set quadrature order appropriate for the distributions and desired accuracy
quadOrder = 15
print(f'mean: {reSurf.getMean(pdfs,quadOrder)}')
v = reSurf.getVariance(pdfs, quadOrder)
print(f'variance: {v[0]}')

# optimize the surrogate
print(f'Optimum: {reSurf.optimize()}')

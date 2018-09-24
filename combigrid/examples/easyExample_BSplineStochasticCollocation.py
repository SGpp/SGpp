import numpy as np
import pysgpp
import time

def f(x):
    return x[0]*x[0]

numDimensions = 1
degree = 5

# set the pointHierarchies to exponentially growing uniform boundary grids in every dimension
pointHierarchies = pysgpp.AbstractPointHierarchyVector()
for d in range(numDimensions):
    pointHierarchies.push_back(pysgpp.CombiHierarchies.expUniformBoundary())
    
# create the BsplinecoefficientGridFunction which calculates the coefficients of the B-spline interpolation of 
# the objective function    
objFct = pysgpp.multiFunc(f)
gf = pysgpp.BSplineCoefficientGridFunction(objFct, pointHierarchies, degree)

# create the CombigridSurrogateModelConfiguration
coefficient_op = pysgpp.CombigridMultiOperation.createExpUniformBoundaryBsplineInterpolation(numDimensions,
                                                                                            objFct,
                                                                                            degree)

config = pysgpp.CombigridSurrogateModelConfiguration()
config.type = pysgpp.CombigridSurrogateModelsType_BSPLINE_STOCHASTIC_COLLOCATION
config.pointHierarchies = pointHierarchies
config.levelManager = pysgpp.RegularLevelManager()
config.degree = degree
config.storage = coefficient_op.getStorage()
config.numDimensions = numDimensions

# create the Bspline Stochasti Collocation from the configuration
bsc = pysgpp.BsplineStochasticCollocation(config)

# create a B-spline interpolation operation, add regular levels and extract the interpolation coefficients and the level structure
operation = pysgpp.CombigridOperation_createExpUniformBoundaryBsplineInterpolation(numDimensions, objFct, degree)
q = 8
operation.getLevelManager().addRegularLevels(q)


config.coefficientStorage = operation.getStorage()
config.levelStructure = operation.getLevelManager().getLevelStructure() 
#print("level structure: ")
#operation.getLevelManager().printLevelStructure(config.levelStructure)
bsc.updateConfig(config) 

# calculate the variance for the first time. This includes the transformation to the hierarchical basis
start = time.time()
variance = bsc.variance()
variancetime = time.time() - start
print("#gp before coarsen {}".format(bsc.numHierarchicalGridPoints()))
print("Variance: {}, Time (trafo+calc):         {}".format(variance, variancetime))

# calculate variance again. This time the transformation has already been done. run time should be much smaller
start = time.time()
variance = bsc.variance()
variancetime = time.time() - start
print("Variance: {}, Time (calc):               {}".format(variance, variancetime))

# coarsen the grid, removing less important grid points and calculate the variance again.
bsc.coarsen(250,0.2)
print("#gp after coarsen {}".format(bsc.numHierarchicalGridPoints()))
start = time.time()
variance = bsc.variance()
variancetime = time.time() - start
print("Variance: {}, Time(calc on coarse grid): {}".format(variance, variancetime))






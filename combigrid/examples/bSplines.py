from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation,\
    CombigridMultiOperation, CombigridTensorOperation
import pysgpp
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, dblquad
from pysgpp.extensions.datadriven.uq.dists import Uniform
from pysgpp.extensions.datadriven.uq.dists.Beta import Beta
from numpy import square

def objectiveFct(x):
    #return np.arctan(50.0 * (x[0] - .35)) + np.pi / 2.0 + 4.0 * x[1] ** 3 + np.exp(x[0] * x[1] - 1.0)
    return np.arctan(50.0 * (x[0] - .35)) #exact: 0.4589

func = pysgpp.multiFunc(lambda x: objectiveFct(x))
numDims = 1
degree = 3
minLevel = 0
maxLevel = 8



grids = pysgpp.AbstractPointHierarchyVector()
evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
for d in range(0,numDims):
    grids.push_back(pysgpp.CombiHierarchies.expUniformBoundary())
    evaluators.push_back(pysgpp.CombiEvaluators.BSplineInterpolation(degree))

gf =pysgpp.BSplineCoefficientGridFunction(func, grids, degree)
exploitNesting = False
storage = pysgpp.CombigridTreeStorage(grids,exploitNesting)
fullGridEval = pysgpp.ScalarFullGridGridBasedEvaluator(storage, evaluators, grids, gf)

def BsplineInterpolation(x,levels):
    params = pysgpp.FloatScalarVectorVector()
    params.push_back(pysgpp.FloatScalarVector(x))
    fullGridEval.setParameters(params)
    result = fullGridEval.eval(levels)
    return result.getValue()
    
for level in range(minLevel,maxLevel+1): 
    levels = [level]
    integral = quad(lambda x: BsplineInterpolation(x,levels), 0,1)
        #meanSquare = dblquad(lambda x,y: BsplineInterpolation(x,y,levels)**2, 0,1,lambda x:0,lambda x:1)
        #variance = meanSquare[0] - mean[0]**2
  
    print"%i  |  integral %f" %(level,integral[0])  

#===============================================================================
#def BsplineInterpolation(x,y, levels):
#    params = pysgpp.FloatScalarVectorVector()
#    params.push_back(pysgpp.FloatScalarVector(x))
#    params.push_back(pysgpp.FloatScalarVector(y))
#    fullGridEval.setParameters(params)
#    result = fullGridEval.eval(levels)
#    return result.getValue()
# for level1 in range(minLevel,maxLevel+1): 
#     for level2 in range(minLevel,maxLevel+1 - level1):
#         levels = [level1, level2]  
#         integral = dblquad(lambda x,y: BsplineInterpolation(x,y,levels), 0,1,lambda x:0,lambda x:1)
#         #meanSquare = dblquad(lambda x,y: BsplineInterpolation(x,y,levels)**2, 0,1,lambda x:0,lambda x:1)
#         #variance = meanSquare[0] - mean[0]**2
#   
#         print"%i %i  |  integral %f" %(level1,level2,integral[0])  
#===============================================================================



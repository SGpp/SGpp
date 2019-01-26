#%matplotlib inline
import matplotlib.pyplot as plt
from sys import path
import numpy as np

#from numpy import *
import scipy as sp
from spatiallyAdaptiveExtendSplit import *
from spatiallyAdaptiveSplit import *
from spatiallyAdaptiveSingleDimension import *
from spatiallyAdaptiveCell import *

from PerformTestCase import *
from Function import *
from ErrorCalculator import *
import math
dim = 2
a = np.zeros(dim)
b = np.ones(dim)
midpoint = np.ones(dim) * 0.5
coefficients = np.array([ 10**0 * (d+1) for d in range(dim)])
f = GenzDiscontinious(border=midpoint,coeffs=coefficients)
f.plot(np.ones(dim)*a,np.ones(dim)*b)
errorOperator=ErrorCalculatorSurplusCell()
errorOperator2=ErrorCalculatorAnalytic()

grid=TrapezoidalGrid(a, b)
adaptiveCombiInstanceSingleDim = SpatiallyAdaptiveSingleDimensions(a, b,grid)
adaptiveCombiInstanceFixed = SpatiallyAdaptiveFixedScheme(a, b,grid)
adaptiveCombiInstanceExtend = SpatiallyAdaptiveExtendScheme(a, b,2,grid,version=0)
adaptiveCombiInstanceCell = SpatiallyAdaptiveCellScheme(a, b,grid)

adaptiveCombiInstanceExtend.performSpatiallyAdaptiv(1,2,f,errorOperator2,10**-2, do_plot=True)

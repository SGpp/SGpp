
# coding: utf-8

# # The combigrid module
# This tutorial contains examples with increasing complexity to introduce you to the combigrid
# module. The combigrid module is quite separated from the other modules. It only refers to the
# base module for things like DataVector and DataMatrix.
# 

# In[1]:


from __future__ import print_function

from itertools import product, combinations, permutations, combinations_with_replacement
#from pysgpp.extensions.datadriven.uq.dists import J, Beta, Uniform
from pysgpp.extensions.datadriven.uq.plot.colors import initialize_plotting_style, load_color, load_font_properties, savefig
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d
import math
import pysgpp

from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib.pyplot as plt
import numpy as np


# ## Test functions
# We define two arbitrary-dimensional functions and choose $d=3$ for later use

# In[2]:


def f(x):
    product = 1.0
    for i in range(x.getSize()):
        product *= math.exp(-x[i])
    return product


def g(x):
    return np.prod([4 * xi * (1 - xi) for xi in x.array()])

d=3


# We have to wrap f in a pysgpp.MultiFunction object to be able to pass a Python function to C++
# and choose $d=3$

# In[3]:


func = pysgpp.multiFunc(f)


# ## Example: Leja quadrature with linear growth of grid points
# Linear growth of the size of the grid wrt. the level

# In[4]:


growthFactor = 2


# Create the operation object that handles the evaluation. The evaluation mode is quadrature,
# so it will approximate the integral of f over $[0, 1]^d$. It uses Leja points with $1 + 2*l$
# points in level l. The level starts from zero, higher level means finer grid.
# Slower growth of the number of points per level means that the total number of points used can
# be controlled better.
# 

# In[5]:


operation = pysgpp.CombigridOperation.createLinearLejaQuadrature(d, func, growthFactor)


# Compute quadrature for level 2, i.e., $|\vec{l}|_1\leq 2$

# In[6]:


result = operation.evaluate(2)
print("Quadrature result: " + str(result) + ", analytical solution: " + str(math.pow(1.0 - 1.0 / math.e, d)))
print("Number of function evaluations: " + str(operation.numGridPoints()))


# ## Example: Polynomial interpolation on nested Clenshaw Curtis grids
# Example is based on interpolation and Clenshaw Curtis grids (nested, exp. growth)

# In[7]:


operation = pysgpp.CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(d, func)


# Evaluate for level 3

# In[8]:


level = 3
evaluationPoint = pysgpp.DataVector([0.1572, 0.6627, 0.2378])
result = operation.evaluate(level, evaluationPoint)
print("Interpolation result: " + str(result) + ", exact function value: " + str(func(evaluationPoint)))


# ### more elaborate: non-regular level choice
# New evaluation point, reset intermediate values computed so far

# In[9]:


evaluationPoint = pysgpp.DataVector([0.4444, 0.6627, 0.2378])
print("Target function value: " + str(func(evaluationPoint)))
operation.setParameters(evaluationPoint)


# The level manager gives access to the active level set

# In[10]:


levelManager = operation.getLevelManager()
levelManager.addRegularLevels(3)


# Evaluate at the specified point

# In[11]:


print("Regular result 1: " + str(operation.getResult()))
print("Total function evaluations: " + str(operation.numGridPoints()))


# Increase regular level to contain at most 50 extra function evaluations.
# 
# Evaluate f automatically and in parallel (4 processes).

# In[12]:


maxNumEvals = 50
numThreads = 4
levelManager.addRegularLevelsByNumPointsParallel(maxNumEvals, numThreads)
print("Regular result 2: " + str(operation.getResult()))
print("Total function evaluations: " + str(operation.numGridPoints()))


# Now adapt with a refinement strategy based on the evaluation point we specified

# In[13]:


operation.setLevelManager(pysgpp.AveragingLevelManager())
levelManager = operation.getLevelManager()
maxNumEvals = 60
levelManager.addLevelsAdaptive(maxNumEvals)
print("Adaptive result: " + str(operation.getResult()))
print("Total function evaluations: " + str(operation.numGridPoints()))


# ### ... and plot the grid

# In[14]:


grid = levelManager.getGridPointMatrix()
gridList = [[grid.get(r, c) for c in range(grid.getNcols())] for r in range(grid.getNrows())]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.scatter(gridList[0], gridList[1], gridList[2], c='r', marker='o')


# ## Example: different grids/operations in different dimensions
# 1. dimension 0: (nested) Chebyshev points ($n = 3^l$ points at level $l$
# 2. dimension 1: Leja points with linear growth $n = 1 + 3l$
# 3. dimension 2: equidistand points with boundary, $n = 2^l+1$

# In[15]:


grids = pysgpp.AbstractPointHierarchyVector()
grids.push_back(pysgpp.CombiHierarchies.expChebyshev())
grids.push_back(pysgpp.CombiHierarchies.linearLeja(3))
grids.push_back(pysgpp.CombiHierarchies.expUniformBoundary())


# Different linear operations per dimension
# 1. dimension 0: polynomial interpolation
# 2. dimension 1: quadrature
# 3. dimension 2: linear interpolation

# In[16]:


evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
evaluators.push_back(pysgpp.CombiEvaluators.polynomialInterpolation())
evaluators.push_back(pysgpp.CombiEvaluators.quadrature())
evaluators.push_back(pysgpp.CombiEvaluators.linearInterpolation())


# now evaluate

# In[17]:


levelManager = pysgpp.WeightedRatioLevelManager()
operation = pysgpp.CombigridOperation(grids, evaluators, levelManager, func)
parameters = pysgpp.DataVector([0.777, 0.14159])
result = operation.evaluate(2, parameters)
print("Result: " + str(result))


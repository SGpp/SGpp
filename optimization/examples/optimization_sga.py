
# coding: utf-8

# # Optimization example

# We look at an example application of the sgpp::optimization module.

# The example interpolates a bivariate test function like the \ref example_tutorial_cpp example.
# However, we use B-splines here instead to obtain a smoother interpolant.
# The resulting sparse grid function is then minimized with the method of steepest descent.
# For comparison, we also minimize the objective function with Nelder-Mead's method.

# First, we import pysgpp and the required modules.

# In[1]:


from __future__ import print_function
import pysgpp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import sys
import numpy as np


# The function $f\colon [0, 1]^d \to \mathbb{R}$ to be minimized
# is called <i>objective function</i> and has to derive from pysgpp.OptScalarFunction.
# In the constructor, we give the dimensionality of the domain
# (in this case $d = 2$).
# The eval method evaluates the objective function and returns the function
# value $f(\vec{x})$ for a given point $\vec{x} \in [0, 1]^d$.

# In[2]:


class ExampleFunction(pysgpp.OptScalarFunction):
    """Example objective function."""
    def __init__(self):
        super(ExampleFunction, self).__init__(2)

    def eval(self, x):
        """Evaluates the function."""
        return np.sin(8.0 * x[0]) + np.sin(7.0 * x[1])


# In[3]:


XX = np.linspace(0,1,17)
(X, Y) = np.meshgrid(XX, XX)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
def ff(x1, x2):
    return np.sin(8.0*x1) + np.sin(7.0*x2)
ax.plot_surface(X, Y, ff(X, Y), rstride=1, cstride=1, cmap='coolwarm')


# We have to disable OpenMP within pysgpp since it interferes with SWIG's director feature.

# In[4]:


pysgpp.omp_set_num_threads(1)

# increase verbosity of the output
pysgpp.OptPrinter.getInstance().setVerbosity(2)


# Here, we set define some parameters: objective function, dimensionality,
# B-spline degree, maximal number of grid points, and adaptivity.

# In[5]:


# objective function
f = ExampleFunction()
# dimension of domain
d = f.getNumberOfParameters()
# B-spline degree
p = 3
# maximal number of grid points
N = 30
# adaptivity of grid generation
gamma = 0.95


# First, we define a grid with modified B-spline basis functions and
# an iterative grid generator, which can generate the grid adaptively.

# In[6]:


grid = pysgpp.Grid.createModBsplineGrid(d, p)
gridGen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma)


# With the iterative grid generator, we generate adaptively a sparse grid.

# In[7]:


if not gridGen.generate():
    print("Grid generation failed, exiting.")
    sys.exit(1)


# Then, we hierarchize the function values to get hierarchical B-spline
# coefficients of the B-spline sparse grid interpolant
# $\tilde{f}\colon [0, 1]^d \to \mathbb{R}$.

# In[8]:


functionValues = gridGen.getFunctionValues()
coeffs = pysgpp.DataVector(len(functionValues))
hierSLE = pysgpp.OptHierarchisationSLE(grid)
sleSolver = pysgpp.OptAutoSLESolver()

# solve linear system
if not sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs):
    print("Solving failed, exiting.")
    sys.exit(1)


# We define the interpolant $\tilde{f}$ and its gradient
# $\nabla\tilde{f}$ for use with the gradient method (steepest descent).
# Of course, one can also use other optimization algorithms from
# sgpp::optimization::optimizer.

# In[9]:


ft = pysgpp.OptInterpolantScalarFunction(grid, coeffs)
ftGradient = pysgpp.OptInterpolantScalarFunctionGradient(grid, coeffs)
gradientDescent = pysgpp.OptGradientDescent(ft, ftGradient)
x0 = pysgpp.DataVector(d)


# The gradient method needs a starting point.
# We use a point of our adaptively generated sparse grid as starting point.
# More specifically, we use the point with the smallest
# (most promising) function value and save it in x0.

# In[10]:


gridStorage = gridGen.getGrid().getStorage()

# index of grid point with minimal function value
x0Index = 0
fX0 = functionValues[0]
for i in range(1, len(functionValues)):
    if functionValues[i] < fX0:
        fX0 = functionValues[i]
        x0Index = i

x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));
ftX0 = ft.eval(x0)

print("x0 = {}".format(x0))
print("f(x0) = {:.6g}, ft(x0) = {:.6g}\n".format(fX0, ftX0))


# We apply the gradient method and print the results.

# In[11]:


gradientDescent.setStartingPoint(x0)
gradientDescent.optimize()
xOpt = gradientDescent.getOptimalPoint()
ftXOpt = gradientDescent.getOptimalValue()
fXOpt = f.eval(xOpt)

print("\nxOpt = {}".format(xOpt))
print("f(xOpt) = {:.6g}, ft(xOpt) = {:.6g}\n".format(fXOpt, ftXOpt))


# For comparison, we apply the classical gradient-free Nelder-Mead method
# directly to the objective function $f$.

# In[12]:


nelderMead = pysgpp.OptNelderMead(f, 1000)
nelderMead.optimize()
xOptNM = nelderMead.getOptimalPoint()
fXOptNM = nelderMead.getOptimalValue()
ftXOptNM = ft.eval(xOptNM)

print("\nxOptNM = {}".format(xOptNM))
print("f(xOptNM) = {:.6g}, ft(xOptNM) = {:.6g}\n".format(fXOptNM, ftXOptNM))


# We see that both the gradient-based optimization of the smooth sparse grid
# interpolant and the gradient-free optimization of the objective function
# find reasonable approximations of the minimum, which lies at
# $(3\pi/16, 3\pi/14) \approx (0.58904862, 0.67319843)$.

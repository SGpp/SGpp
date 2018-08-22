
# coding: utf-8

# In[1]:

from __future__ import print_function


# In[2]:


from pysgpp import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math


# # Interpolation
# of a function $f: \Omega:=[0,1]^2\to\mathbb{R}$ with $\delta\Omega=0$

# In[3]:


def f(x1, x2):
    return 16*x1*(1-x1**2)*x2*(1-x2**2)
# and plot
XX = np.linspace(0,1,17)
(X, Y) = np.meshgrid(XX, XX)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, f(X, Y), rstride=1, cstride=1, cmap='coolwarm')


# In[4]:


dim = 2
grid = Grid.createLinearGrid(dim)
#grid = Grid.createPolyGrid(dim, 3)
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)


# Separation of concerns:
# GridStorage object allows to provide different implementations

# In[5]:


gridStorage = grid.getStorage()
print("Dim    =", gridStorage.getDimension())
print("|grid| =", gridStorage.getSize())


# Set vector of function values

# In[6]:


alpha = DataVector(gridStorage.getSize())
for i in range(gridStorage.getSize()):
    gp = gridStorage.getPoint(i)
    alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))
print(alpha)


# ### Hierarchize (inplace)
# 
# Factory method chooses best available implementation for given grid

# In[7]:


opHier = createOperationHierarchisation(grid)
opHier.doHierarchisation(alpha)
print(alpha)


# ## plotting

# In[8]:


coords = np.array([[gridStorage.getPoint(i).getStandardCoordinate(0), gridStorage.getPoint(i).getStandardCoordinate(1)] for i in range(gridStorage.getSize())])
fig = plt.figure()
plt.xlim([0,1])
plt.ylim([0,1])
plt.plot(coords[:,0], coords[:,1], 'b.')


# ## evaluate

# In[9]:


opEval = createOperationEval(grid)
p = DataVector([0.3, 0.3])
opEval.eval(alpha, p)


# ## plot new function

# In[10]:


xy = DataMatrix(list(zip(X.flatten(), Y.flatten())))
z = DataVector(xy.getNrows())
opMEval = createOperationMultipleEval(grid, xy)
opMEval.mult(alpha, z)
res = z.array()
res = res.reshape((len(X), len(X[0])))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, res, rstride=1, cstride=1, cmap='coolwarm')


# ## now adapt

# In[11]:


N_old = gridStorage.getSize()
gridGen.refine(SurplusRefinementFunctor(alpha, 1))
N = gridStorage.getSize()
print(N_old, "->", N)


# In[12]:


coords = np.array([[gridStorage.getPoint(i).getStandardCoordinate(0), gridStorage.getPoint(i).getStandardCoordinate(1)] for i in range(gridStorage.getSize())])
fig = plt.figure()
plt.xlim([0,1])
plt.ylim([0,1])
plt.plot(coords[:,0], coords[:,1], 'b.')


# compute new coefficients

# In[13]:


alpha.resizeZero(N)
for i in range(N_old, N):
    pos = [gp.getStandardCoordinate(0), gp.getStandardCoordinate(1)]
    alpha[i] = f(pos[0], pos[1]) - opEval.eval(alpha, DataVector(pos))


# In[14]:


opMEval.mult(alpha, z)
res = z.array()
res = res.reshape((len(X), len(X[0])))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, res, rstride=1, cstride=1, cmap='coolwarm')


# ## and now quadrature

# In[15]:


opQ = createOperationQuadrature(grid)
opQ.doQuadrature(alpha)

plt.show()


import functions
import numpy as np
from argparse import ArgumentParser
import ipdb
import os
import time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import cPickle as pickle
import functions
import numpy as np
import pysgpp

from functions import objFuncSGpp as objFuncSGpp

# f = functions.getFunction('anugaTime', 3)
# y = f.eval([0.8] * 3)
# print(y)
# f.cleanUp()

########################### create random validation data ###############################
dim = 2
num = 10000
r = np.random.rand(num, dim)
path = '/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA/Values'
np.savetxt(os.path.join(path, 'x{}D.txt'.format(dim)), r, delimiter=',') 
f = functions.getFunction('anuga', dim)
values = {}
for i in range(num):
    print("{}:".format(i))
    x = r[i, :]
    values[tuple(x)] = f.eval(x)
          
with open(os.path.join(path, 'values{}D.pkl'.format(dim)), "wb") as f:
            pickle.dump(values, f)
# don't clean up. Otherwise the precalculations file will blow up with unnecessary random values
# f.cleanUp()

######################### plot exemplary timelines ########################################
# dim = 2
# path = '/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA/Values'
# with open(os.path.join(path, 'values{}.pkl'.format(dim)), 'rb') as f:
#     values = pickle.load(f)
#  
# a = values[(9.024547692881779160e-02, 1.971720849876249515e-01)]
# b = values[ (7.702299560766515674e-01, 5.425699119831074446e-01)]
# c = values[(9.427703859143611309e-01, 6.861858146040415996e-01)]
#  
# plt.plot(range(len(a[0])), a[0])
# plt.plot(range(len(b[0])), b[0])
# plt.plot(range(len(c[0])), c[0])
# plt.show()

######################### surf 2D ANUGA ########################################
# dim = 2
# f = functions.getFunction('anuga', dim)
# npts = 10
# x = np.linspace(0, 1, npts)
# X, Y = np.meshgrid(x, x)
# Z = np.zeros(np.shape(X))
# for i in range(npts):
#     for j in range(npts):
#         Z[i, j] = f.eval([X[i, j], Y[i, j]])
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z, cmap=cm.viridis)
# f.cleanUp()
# plt.show()

############################# create and plot timeline surrogate######################################
# model = 'anugaTime'
# dim = 3
# pyFunc = functions.getFunction(model, dim)
# objFunc = objFuncSGpp(pyFunc)
# gridType = "nakbsplineextended"
# degree = 3
# reSurf = pysgpp.SparseGridResponseSurfaceBspline(objFunc,
#                                                  pysgpp.Grid.stringToGridType(gridType),
#                                                  degree)
# level = 4
# reSurf.regular(level) 
# point = pysgpp.DataVector(dim, 0.5)
# 
# x = np.linspace(0, 1, 100)
# F = [0] * len(x)
# for i in range(len(x)):
#     point[2] = x[i]
#     F[i] = reSurf.eval(point)
#     
# plt.plot(x, F)
# plt.show()

from argparse import ArgumentParser
import colorsys
import ipdb
import itertools
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
import os
import scipy

from matplotlib import cm
from matplotlib.pyplot import gca

import active_subspaces as ac
import asFunctions
import cPickle as pickle
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pn
import pysgpp


def MsplineWiki(n, i, x, xi):
#     print("degree: {} index {} x {} xi {}".format(n, l, x, xi))
#     print("{} {} {} {}".format(n, i, xi[i + n], xi[i]))
    if n == 1:
        if xi[i] <= x and x < xi[i + 1]:
            return 1.0 / (xi[i + 1] - xi[i])
        else:
            return 0
    else:
        if (xi[i + n] - xi[i]) != 0:
            return (n * ((x - xi[i]) * MsplineWiki(n - 1, i, x, xi) + (xi[i + n] - x) * MsplineWiki(n - 1, i + 1, x, xi))) / ((n - 1) * (xi[i + n] - xi[i]))
        else:
             return 0.0


def corners(dim):
    corners = np.zeros((dim, 2 ** dim))
    if dim == 1:
        corners[0, 0] = 1
        corners[0, 1] = 0
        return corners
    else:
        jump = 2 ** (dim - 1)
        for i in range(dim):
            j = 0
            while j < 2 ** dim:
                for n in range(jump):
                    if j + n >= 2 ** dim:
                        break
                    corners[i, j + n] = 1
                j = j + 2 * jump
            jump = jump / 2
        return corners

# corresponds to sin(sum(x))
# W1 = ones(dim)/sqrt(dim)
# def func1D(x, dim):
#     return np.sin(dim / np.sqrt(dim) * x)


# corresponds to sin(alpha *sum(x) + 1) / (alpha *sum(x) + 1)
# W1 = ones(dim)/sqrt(dim)
def func1D(x, dim, alpha=1.5):
    t = x * dim * alpha / np.sqrt(dim) + 1
    return np.sin(t) / t


dim = 8
W1 = np.ones(dim) / np.sqrt(dim)  # function dependent!!!
perm = range(dim)
permutations = list(itertools.permutations(perm))
projectedCorners = np.zeros(shape=(dim + 1, int(scipy.misc.factorial(dim))))
for i in range(len(permutations)):
    V = np.zeros((dim, dim + 1))
    for k in range(dim):
        for l in range(k, dim):
            V[permutations[i][l], k] = 1
    projectedCorners[:, i] = W1.transpose().dot(V)
 
leftBound = np.min(projectedCorners)
rightBound = np.max(projectedCorners)

print("bounds: [{}, {}]".format(leftBound, rightBound))
 
 
def Vol(x):  
    vol = 0 
    for i in range(int(scipy.misc.factorial(dim))): 
        xi = projectedCorners[:, i]
        xi.sort()
        vol += MsplineWiki(len(xi) - 1, 0, x, xi)
    return vol


quadArg = lambda  x: func1D(x, dim) * Vol(x)

integral = scipy.integrate.quad(quadArg, leftBound, rightBound)
print(integral[0] / scipy.misc.factorial(dim))

# X = np.linspace(leftBound, rightBound, 100)
# M = np.zeros(len(X))
# F = np.zeros(len(X))
# FV = np.zeros(len(X))
# for i in range(len(X)):
# #     M[i] = Vol(X[i]) 
#     F[i] = func1D(X[i], dim)
# #     FV[i] = func1D(X[i], dim) * Vol(X[i])
# # plt.plot(X, M, 'b')
# plt.plot(X, F, 'g')
# # plt.plot(X, FV, 'yellow')
# plt.show()


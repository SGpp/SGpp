from mpl_toolkits.mplot3d import Axes3D
import os
import sys        
import time

from matplotlib import cm

import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def xpowplus(x, n):
    if x >= 0:
        return x ** (n)
    else:
        return 0.0

    
def w(n, i, v, xi):
    res = 1.0
    for j in range(i, n + i + 1):
        if xi[v] == xi[j]:
            continue
        res *= (xi[v] - xi[j])
    return res

    
# "for multiple knots use the formula for divided differences for multiple arguments"
def Mspline(n, i, x, xi):
    res = 0.0
    for v in range(i, n + i + 1):
        res += (n * xpowplus(xi[v] - x, n - 1)) / w(n, i, v, xi)
    return res


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


def Bspline(n, i, x, xi):
    if n == 0:
        if xi[i] <= x and x < xi[i + 1]:
            return 1.0
        else:
            return 0.0
    else:
        left = 0.0
        right = 0.0
        if (xi[i + n] - xi[i]) != 0:
            left = (x - xi[i]) / (xi[i + n] - xi[i]) * Bspline(n - 1, i, x, xi)
        if (xi[i + n + 1] - xi[i + 1]) != 0:
            right = (xi[i + n + 1] - x) / (xi[i + n + 1] - xi[i + 1]) * Bspline(n - 1, i + 1, x, xi)
        return left + right


xi = [0, 1, 1, 1, 2]
X = np.linspace(xi[0], xi[-1], 100)
M = np.zeros(np.shape(X))
MW = np.zeros(np.shape(X))
B = np.zeros(np.shape(X))
index = 0
degree = 4
for l in range(len(X)):
    M[l] = Mspline(degree, index, X[l], xi)
    MW[l] = MsplineWiki(degree, index, X[l], xi)
    B[l] = Bspline(degree - 1, index, X[l], xi)
      
# print("diff MW {}".format(np.linalg.norm(M - MW)))
# print("diff B {}".format(np.linalg.norm(M - B)))
plt.plot(X, M, 'b')
plt.plot(X, MW, 'g+')
plt.plot(X, B, 'r*')
plt.show()

# MsplineWiki(4, 0, 0.3, xi)


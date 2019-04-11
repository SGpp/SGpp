import pysgpp
import math
import sys
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

import functions


def anugaError(reSurf, objFunc):
    dim = objFunc.getDim()
    if "anugaTime" in objFunc.getName():
        dim -= 1
    # load validation values 
    path = '/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA/Values'
    parameters = np.loadtxt(os.path.join(path, 'x{}D.txt'.format(dim)), delimiter=',') 
    with open(os.path.join(path, 'values{}D.pkl'.format(dim)), 'rb') as f:
        values = pickle.load(f)
    num = np.shape(parameters)[0]
    
    if "anugaTime" in objFunc.getName():
        print("Implement this")  
        return 777      
    else:  # "anuga"
        err = 0
        v = pysgpp.DataVector(dim)
        for i in range(num):
            for d in range(dim):
                v[d] = parameters[i, d]
                y = values[tuple(parameters[i, :])]
            err += (reSurf.eval(v) - np.max(y)) ** 2
        return np.sqrt(err)


pysgpp.omp_set_num_threads(1)
pysgpp.OptPrinter.getInstance().setVerbosity(-1)

dim = 4
pyFunc = functions.getFunction('anuga', dim)
objFunc = functions.objFuncSGppSign(pyFunc)

degree = 3
nPoints = 20
gridType = "nakbsplineextended"

reSurf = pysgpp.SparseGridResponseSurfaceBspline(objFunc,
                                                 pysgpp.Grid.stringToGridType(gridType),
                                                 degree)
gamma = 0.95
reSurf.ritterNovak(nPoints, gamma)

grid = reSurf.getGrid()
coeffs = reSurf.getCoefficients()
ft = pysgpp.OptInterpolantScalarFunction(grid, coeffs)
ftGradient = pysgpp.OptInterpolantScalarFunctionGradient(grid, coeffs)

xOpt = reSurf.optimize()
# we maximized, so fix the sign again
fXOpt = -objFunc.eval(xOpt)
approxFXOpt = -reSurf.eval(xOpt)

print("l2 error: {}".format(anugaError(reSurf, objFunc)))

print("\nxOpt = {}".format(xOpt))
print("f(xOpt) = {:.6g} approxF(xOpt) = {:.6g}\n".format(fXOpt, approxFXOpt))

pdfs = pysgpp.DistributionsVector(dim, pysgpp.DistributionUniform(0, 1))
quadOrder = 10
# we maximized, so fix the sign again
mean = -reSurf.getMean(pdfs, quadOrder)
var = reSurf.getVariance(pdfs, quadOrder)

print("mean: {}".format(mean))
print("var: {}\n".format(var))

pyFunc.cleanUp()

# plot 2D grid
# gridStorage = grid.getStorage()
# xGrid = [0] * gridStorage.getSize()
# yGrid = [0] * gridStorage.getSize()
# for i in range(gridStorage.getSize()):
#     p = gridStorage.getPointCoordinates(i)
#     xGrid[i] = p[0]
#     yGrid[i] = p[1]
# plt.scatter(xGrid, yGrid)
# plt.show()

####################### Nelder-Mead for comparison ########################
# # For comparison, we apply the classical gradient-free Nelder-Mead method
# # directly to the objective function \f$f\f$.
# printLine()
# print "Optimizing objective function (for comparison)...\n"
# nelderMead = pysgpp.OptNelderMead(objFunc, 1000)
# nelderMead.optimize()
# xOptNM = nelderMead.getOptimalPoint()
# fXOptNM = nelderMead.getOptimalValue()
# ftXOptNM = ft.eval(xOptNM)
# 
# print "\nxOptNM = {}".format(xOptNM)
# print "f(xOptNM) = {:.6g}, ft(xOptNM) = {:.6g}\n".format(fXOptNM, ftXOptNM)
# 
# printLine()
# print "\nsgpp::optimization example program terminated."

# # The example program outputs the following results:
# # \verbinclude optimization.output.txt
# #
# # We see that both the gradient-based optimization of the smooth sparse grid
# # interpolant and the gradient-free optimization of the objective function
# # find reasonable approximations of the minimum, which lies at
# # \f$(3\pi/16, 3\pi/14) \approx (0.58904862, 0.67319843)\f$.

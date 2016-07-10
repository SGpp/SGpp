#!/usr/bin/python

import pysgpp
import math
import sys
import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from matplotlib2tikz import save as tikz_save


class ExampleFunction(pysgpp.OptScalarFunction):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(ExampleFunction, self).__init__(2)

    def eval(self, x):
        # return 0.1*(x*15-8)**2 - 5*np.cos(x*15-8)
        # return target_func(x)
        # return shcb(x)
        # return easom(x)
        # return michalewicz(x)
        return schwefel(x)

def target_func(x):
    grid = pysgpp.Grid.createBsplineGrid(1, 3)
    gridGen = grid.getGenerator().regular(1)
    gridStorage = grid.getStorage()
    alpha = pysgpp.DataVector(gridStorage.getSize())
    hierarch = pysgpp.createOperationMultipleHierarchisation(grid)
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        # x = (gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))
        y = gp.getStandardCoordinate(0)
        # alpha[i] = oscill_genz(d, x)
        # alpha[i] = (0.1*(x[0]*15-8)**2 - 5*np.cos(x[0]*15-8))*2*y
        # alpha[i] = michalewicz(pysgpp.DataVector(x))*2*y
        alpha[i] = easom(pysgpp.DataVector(x))*2*y
        # alpha[i] = 1
    hierarch.doHierarchisation(alpha)
    
    # Quadratur mit f(x) multiplizieren 
    # for i in range(gridStorage.getSize()):
        # gp = gridStorage.getPoint(0)
        # alpha[i] = alpha[i] * f(gp.getCord(0))
    return pysgpp.createOperationQuadrature(grid).doQuadrature(alpha)
    
def shcb(x):
    return x[0]**2 * (4 - 2.1*x[0]**2 + (x[0]**4)/3) + x[0]*x[1] + 4*x[1]**2*(x[1]**2 - 1) 
def oscill_genz(d, x):
    u = 0.5
    a = [5, 5]
    return np.cos(2*np.pi*u+sum([a[i]*x[i] for i in range(d)]))

def michalewicz(x):
    return -np.sin(5*x[0])*np.sin((x[0]*5)**2/np.pi)**20-np.sin(5*x[1])*np.sin(2*(x[1]*5)**2/np.pi)**20
          
def easom(x):
    x[0] = (x[0]+100)/200
    x[1] = (x[1]+100)/200
    return -np.cos(x[0])*np.cos(x[1])*np.exp(-(x[0]-np.pi)**2-(x[1]-np.pi)**2)

def schwefel(x):
    for x_i in x:
        x_i = (x_i + 500) / 1000
    return -sum([x_t*np.sin(abs(x_t)**0.5) for x_t in x] )
        
def grid2d(start, end, num=50):
    """Create an 2D array where each row is a 2D coordinate.
    np.meshgrid is pretty annoying!
    """
    dom = np.linspace(start, end, num)
    X0, X1 = np.meshgrid(dom, dom)
    return np.column_stack([X0.flatten(), X1.flatten()])

def gridToName(i):
    types = {0 : "Linear",
             1 : "LinearStretched",
             2 : "LinearL0Boundary",
             3 : "LinearBoundary",
             4 : "LinearStretchedBoundary",
             5 : "LinearTruncatedBoundary",
             6 : "ModLinear",
             7 : "Poly (Grad 3)",
             8 : "PolyBoundary",
             9 : "ModPoly",
             10 : "ModWavelet",
             11 : "ModBspline (Grad 3)",
             12 : "Prewavelet",
             13 : "SquareRoot",
             14 : "Periodic",
             15 : "LinearClenshawCurtis",
             16 : "B-Spline (Grad 3)",
             17 : "B-SplineBoundary (Grad 3)",
             18 : "B-SplineClenshawCurtis (Grad 3)",
             19 : "Wavelet",
             20 : "WaveletBoundary",
             21 : "FundamentalSpline",
             22 : "ModFundamentalSpline",
             23 : "ModBsplineClenshawCurtis",
             24 : "LinearStencil",
             25 : "ModLinearStencil"}
    return types[i]

def printLine():
    print "----------------------------------------" + \
          "----------------------------------------"


def optimize():
    # disable multi-threading
    pysgpp.omp_set_num_threads(1) # increase output verbosity
    pysgpp.OptPrinter.getInstance().setVerbosity(2)
    print "sgpp::optimization example program started.\n"
    
    # objective function
    f = ExampleFunction()
    # dimension of domain
    d = f.getNumberOfParameters()
    # B-spline degree
    p = 3
    # maximal number of grid points
    # adaptivity of grid generation
    gamma = 0.95 
    grids = (pysgpp.Grid.createBsplineGrid(d, 3), pysgpp.Grid.createBsplineBoundaryGrid(d, 3), pysgpp.Grid.createModBsplineGrid(d, 3))# pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3))
    N_vals =  range(100, 3000, 300)
    for grid in grids:
        
        printLine()
        print gridToName(grid.getType())
        printLine()
        errors = []
        for N in N_vals: 
            grid.getStorage().clear()
            gridGen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma)

            # #############################################################################
            # GRID GENERATION
            # #############################################################################

            printLine()
            print "Generating grid...\n"


            if not gridGen.generate():
                print "Grid generation failed, exiting."
                sys.exit(1)

            # #############################################################################
            # HIERARCHIZATION
            # #############################################################################

            printLine()
            print "Hierarchizing...\n"
            functionValues = gridGen.getFunctionValues()
            coeffs = pysgpp.DataVector(len(functionValues))
            hierSLE = pysgpp.OptHierarchisationSLE(grid)
            sleSolver = pysgpp.OptAutoSLESolver()

            # solve linear system
            if not sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs):
                print "Solving failed, exiting."
                sys.exit(1)

            # #############################################################################
            # OPTIMIZATION OF THE SMOOTH INTERPOLANT
            # #############################################################################

            printLine()
            print "Optimizing smooth interpolant...\n"
            ft = pysgpp.OptInterpolantScalarFunction(grid, coeffs)
            ftGradient = pysgpp.OptInterpolantScalarFunctionGradient(grid, coeffs)
            gradientMethod = pysgpp.OptGradientDescent(ft, ftGradient)
            x0 = pysgpp.DataVector(d)

            # determine best grid point as starting point
            gridStorage = gridGen.getGrid().getStorage()

            # index of grid point with minimal function value
            x0Index = 0
            fX0 = functionValues[0]
            for i in range(1, len(functionValues)):
                if functionValues[i] < fX0:
                    fX0 = functionValues[i]
                    x0Index = i

            for t in range(d):
                x0[t] = gridStorage.getPoint(x0Index).getStandardCoordinate(t)

            ftX0 = ft.eval(x0)

            print "x0 = {}".format(x0)
            print "f(x0) = {:.6g}, ft(x0) = {:.6g}\n".format(fX0, ftX0)

            gradientMethod.setStartingPoint(x0)
            gradientMethod.optimize()
            xOpt = gradientMethod.getOptimalPoint()
            ftXOpt = gradientMethod.getOptimalValue()
            fXOpt = f.eval(xOpt)
            # errors.append(fXOpt - (-1.8013)) #michalewicz
            # errors.append(fXOpt - (-1.031628)) #SHCB
            errors.append(fXOpt - (-1)) #Easom
            
            print "\nxOpt = {}".format(xOpt)
            print "f(xOpt) = {:.6g}, ft(xOpt) = {:.6g}\n".format(fXOpt, ftXOpt)

        plt.plot(N_vals, errors, label=gridToName(grid.getType()))
        # plt.loglog(N_vals, errors, basex=2, basey=2, label=gridToName(grid.getType()))
    plt.xlabel('Grid Points')
    plt.ylabel('Error')
    plt.legend(loc=1)
    plt.title("Optimization-Error SHCB-Funktion") 
    plt.show()
   
def genz_error(sol):
    a = (5,5)
    exact_sol = (np.cos(a[0] + a[1]) - np.cos(a[0]) - np.cos(a[1])+1)/(a[0]*a[1])
    return abs(exact_sol-sol)

def l2_error(d, f, f_sg):
    n = 10000
    s = 0
    for i in range(n):
        x = pysgpp.DataVector(d)
        for j in range(d):
            x[j] = random.uniform(0,1)
        s += (f_sg(x)-f(x))**2
    return (s/n)**0.5
    
def integrate():
    d = 2
    grids = (pysgpp.Grid.createLinearGrid(d), pysgpp.Grid.createPolyGrid(d,3), pysgpp.Grid.createBsplineGrid(d, 3), pysgpp.Grid.createBsplineBoundaryGrid(d, 3), pysgpp.Grid.createModBsplineGrid(d, 3), pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3))
    for grid in grids:
        printLine()
        print gridToName(grid.getType())
        printLine()
        opQ = pysgpp.createOperationQuadrature(grid)
        if grid.getType() == 7:
            hierarch = pysgpp.createOperationHierarchisation(grid)
        else:
            hierarch = pysgpp.createOperationMultipleHierarchisation(grid)
        # system = pysgpp.OptHierarchisationSLE(grid)
        # solver = pysgpp.OptUMFPACK()
        num_gridPoints = []
        errors = []
        for l in range(1,9):
            grid.getStorage().clear()
            gridGen = grid.getGenerator().regular(l)
            gridStorage = grid.getStorage()
            alpha = pysgpp.DataVector(gridStorage.getSize())
            opeval = pysgpp.createOperationNaiveEval(grid)
            for i in range(gridStorage.getSize()):
                gp = gridStorage.getPoint(i)
                x = (gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))
                alpha[i] = oscill_genz(d,x)
                # alpha[i] = target_func(x)
                # alpha[i] = michalewicz(x)
                # alpha[i] = 1
                # print(alpha[i])
            hierarch.doHierarchisation(alpha)
            # b = pysgpp.DataVector(alpha)
            # solver.solve(system, b, alpha)
  
            f_sg = lambda point: opeval.eval(alpha, point)
            # errors.append(l2_error(d, lambda point: oscill_genz(d, point), f_sg) )
            # errors.append(l2_error(d, michalewicz, f_sg) )
            sol = opQ.doQuadrature(alpha)
            errors.append(genz_error(sol))
            # print errors[-1]
            num_gridPoints.append(gridStorage.getSize())
            
        
        plt.loglog(num_gridPoints, errors, basex=10, basey=10, label=gridToName(grid.getType()))
        # plt.plot(range(1,9), errors, label=gridToName(grid.getType()))
    plt.xlabel('#Gitterpunkte')
    plt.ylabel('Fehler')
    plt.legend(loc=3)
    plt.title("Quadratur-Fehler Genz-Funktion") 
    plt.show()
    # tikz_save('mytikz.tex');

def test():
    pysgpp.OptPrinter.getInstance().setVerbosity(10)
    grid = pysgpp.Grid.createBsplineClenshawCurtisGrid(1, 3)
    gridGen = grid.getGenerator().regular(6)
    gridStorage = grid.getStorage()
    alpha = pysgpp.DataVector(gridStorage.getSize())
    system = pysgpp.OptHierarchisationSLE(grid)
    solver = pysgpp.OptArmadillo()
    for i in range(gridStorage.getSize()):
        alpha[i] = 1
    b = pysgpp.DataVector(alpha)
    
    n = system.getDimension()
    print "[" + "\n".join(["[" + ", ".join([str(system.getMatrixEntry(i, j))
                                        for j in range(n)]) + "]"
                       for i in range(n)]) + "]"
    solver.solve(system, b, alpha)
    
def plot_genz():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X = np.arange(0, 1, 0.01)
    Y = np.arange(0, 1, 0.01)
    X, Y = np.meshgrid(X, Y)
    Z = oscill_genz(2, [X,Y])
    ax.plot_surface(X, Y, Z, rstride=4, cstride=4, cmap=cm.coolwarm)
    plt.show()
# test()
# plot_genz()
integrate()
# optimize()
# X_grid = np.outer(np.linspace(0, 1, 50), np.ones(50))
# Y_grid = X_grid.copy().T
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# col_stack = np.column_stack((X_grid,Y_grid))
# X_grid = grid2d(0, 1, 50)
# print X_grid[:, 1]
# Z_grid = np.zeros(2500)
# for i, (x,y) in enumerate(np.column_stack((X_grid,Y_grid))):
# for i,x in enumerate(X_grid):
    # Z_grid[i] = shcb([x,y]) 
    # Z_grid[i] = target_func(x) 
    # Z_grid[i] =  1
# Z_grid = np.log(shcb([X_grid,Y_grid])) + 2
# Z_grid = target_func((X_grid, Y_grid))
# Z_grid = michalewicz((X_grid, Y_grid))
# ax.plot_surface(X_grid[:, 0], X_grid[:, 1], y_grid, rstride=4, cstride=4, color='b')

# ax.plot_wireframe(X_grid[:, 0], X_grid[:, 1], Z_grid)
# xs = np.linspace(0,1,1000)
# plt.plot(xs, [target_func([x]) for x in xs])
# plt.show()


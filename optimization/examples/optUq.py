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
from pysgpp.extensions.datadriven.uq.plot import plot2d
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import checkInterpolation

class ExampleFunction(pysgpp.OptScalarFunction):
  int_grid_level = 3
  def __init__(self):
    super(ExampleFunction, self).__init__(3)
    
  def eval(self, x):
    return self.target_func(x)
    # return beale(x)
    # return shcb(x)
    # return easom(x)
    # return michalewicz(x)   # 
    # return griewank(x)
    # return schwefel(x)
    # return sin_spitze(x)

  def target_func(self,x):
    grid = pysgpp.Grid.createBsplineClenshawCurtisGrid(1, 3)
    gridGen = grid.getGenerator().regular(self.int_grid_level)
    gridStorage = grid.getStorage()
    alpha = pysgpp.DataVector(gridStorage.getSize())
    hierarch = pysgpp.createOperationMultipleHierarchisation(grid)
    for i in range(gridStorage.getSize()):
      gp = gridStorage.getPoint(i)
      y = gridStorage.getCoordinates(gp)
      # alpha[i] = oscill_genz(d, x)
      # alpha[i] = (0.1*(x[0]*15-8)**2 - 5*np.cos(x[0]*15-8))*2*y
      # alpha[i] = michalewicz(pysgpp.DataVector(x))*2*y #
      # alpha[i] = easom(pysgpp.DataVector(x))*2*y
      # alpha[i] = 1
      # alpha[i] = schwefel([x[0], x[1], x[2], y[0]])
      # alpha[i] = revenue(x[0], y[0])
      # alpha[i] = beale([x[0], y[0]])
      # alpha[i] = eggholder([x[0], y[0]])
      # alpha[i] = griewank([x[0], x[1], x[2], y[0]])
      alpha[i] = sin_spitze([x[0], x[1], x[2], y[0]])
    hierarch.doHierarchisation(alpha)
    # Quadratur mit f(x) multiplizieren 
    return pysgpp.createOperationQuadrature(grid).doQuadrature(alpha) 
    
def revenue(z, D):
  s = 30.0
  v = 30.0
  my = 60.0
  scaling = 2*v
  r = 4.0
  c = 2.0
  z *= scaling
  a = s + z 
  D = (D * scaling) + (my-v)
  assert s <= a
  return (r - c) * a + c*s - (r - 0.5*c) * max(a - D, 0) 

def expected_revenue(z):
  s = 30.0
  v = 30.0
  my = 60.0
  scaling = 2*v
  r = 4.0
  c = 2.0
  z *= scaling
  a = s + z
  assert s <= a
  expected = 0
  if (my - v < a and a <= my + v):
    expected = ((a - my + v)**2)/(4*v)
  elif(a > my + v):
    expected = a - my
  # print my + (v*(2*r - 3*c))/(2*r - c) # opt a
  print (my + (v*(2*r - 3*c))/(2*r - c) - s) / scaling 
  return (r - c) * a + c*s - (r - 0.5*c) * expected
   
def sin_kuppel(x):
  return np.sin(np.pi*x[0])*np.sin(np.pi*x[1])
    
def gamma(n):
  res = 1
  for i in range(1, n):
    res *= i
  return float(res)
    
def beta(p, q):
  return (gamma(p)*gamma(q)) / gamma(p+q)

def f_beta(x):
  p = 3
  q = 4
  return (1.0/beta(p, q)) * x**(p-1) * (1-x)**(q-1) 

def shcb(x):
  x[0] = x[0]*10 - 5
  x[1] = x[1]*10 - 5
  return x[0]**2 * (4 - 2.1*x[0]**2 + (x[0]**4)/3) + x[0]*x[1] + 4*x[1]**2*(x[1]**2 - 1) 
    
def oscill_genz(d, x):
  u = 0.5
  a = [2, 2]
  return np.cos(2*np.pi*u+sum([a[i]*x[i] for i in range(d)]))

def michalewicz(x):
  x[0] = x[0]*5
  x[1] = x[1]*5
  return -np.sin(x[0])*np.sin(x[0]**2/np.pi)**20-np.sin(x[1])*np.sin(2*x[1]**2/np.pi)**20

def eggholder(x):
  x[0] = x[0]*1024 - 512
  x[0] = x[1]*1024 - 512
  return -(x[1]+47)*np.sin(abs(x[0]/2 + x[1] + 47)**0.5)-x[0]*np.sin(abs(x[0]-x[1]-47)**0.5)
def eggcrate(x):
  x[0] = x[0]*10 - 5
  x[1] = x[1]*10 - 5
  return x[0]**2 + x[1]**2 + 25*(np.sin(x[0])**2 + np.sin(x[1])**2)
          
def easom(x):
  x[0] = x[0]*200 - 100
  x[1] = x[1]*200 - 100
    
  return -np.cos(x[0])*np.cos(x[1])*np.exp(-(x[0]-np.pi)**2-(x[1]-np.pi)**2)
  
def sin_spitze(x):
  d = len(x)
  prod = 1.0
  for i in range(d):
    prod *= x[i]*np.sin(5*np.pi*x[i]+2*np.pi*x[d-i-1])
  return -prod

def griewank(x):  #falsch
  #np.linalg.norm is buggy for pysgpp.DataVector
  d = len(x)
  prod = 1.0
  norm = 0.0
  for i in range(d):
    x[i] = x[i]*1200.0 - 601.0
    prod *= np.cos(x[i] / (i+1.0)**0.5 ) 
    norm += x[i]**2
  norm = norm**0.5
  return 1.0 + norm**2 / 4000.0 - prod
    
def beale(x):
  x[0] = x[0] * 10 - 5
  x[1] = x[1] * 10 - 5
  return (1.5 - x[0]*(1-x[1]))**2+(2.25 - x[0]*(1-x[1]**2))**2+(2.625-x[0]*(1-x[1]**3))**2

  
def schwefel(x):
    for i in range(len(x)):
        x[i] = x[i] * 1000 - 500
    return -sum([x[i]*np.sin(abs(x[i])**0.5) for i in range(len(x))] )
        
def grid2d(start, end, num=50):
    """Create an 2D array where each row is a 2D coordinate.
    np.meshgrid is pretty annoying!
    """
    dom = np.linspace(start, end, num)
    X0, X1 = np.meshgrid(dom, dom)
    return np.column_stack([X0.flatten(), X1.flatten()])

def gridToName(i, p):
    types = {0 : "Linear",
             1 : "LinearStretched",
             2 : "LinearL0Boundary",
             3 : "LinearBoundary (Rand)",
             4 : "LinearStretchedBoundary", 5 : "LinearTruncatedBoundary",
             6 : "ModLinear",
             7 : "Poly (Grad " + str(p) + ")",
             8 : "PolyBoundary (Grad 3) (Rand)",
             9 : "ModPoly",
             10 : "ModWavelet",
             11 : "ModBspline (Grad 3) (Rand)",
             12 : "Prewavelet",
             13 : "SquareRoot",
             14 : "Periodic",
             15 : "LinearClenshawCurtis",
             16 : "B-Spline (Grad 3)",
             17 : "B-SplineBoundary (Grad " + str(p) + ") (Rand)",
             18 : "B-SplineClenshawCurtis (Grad " + str(p) + ") (Rand)",
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
    N_vals =  range(100, 3301, 800)
    grids = [pysgpp.Grid.createBsplineGrid(d, 3),
             pysgpp.Grid.createBsplineBoundaryGrid(d, 3),
             # pysgpp.Grid.createModBsplineGrid(d,p),
             # pysgpp.Grid.createModBsplineClenshawCurtisGrid(d, 3),
             pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3)]
    
    fig=plt.figure()
    ax=fig.add_subplot(111)
    for grid in grids:
        printLine()
        print gridToName(grid.getType(), p)
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

            x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));
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
            errors.append(fXOpt) #Griewank, Beale
            # errors.append(fXOpt - -160) # revenue
            # errors.append(fXOpt - (-837.966) ) # SChwefel 3d opt
            # errors.append(fXOpt - (-1)) #Easom
            # errors.append(fXOpt - d*(-418.9829)) #Schwefel
            
            print "\nxOpt = {}".format(xOpt)
            print "f(xOpt) = {:.6g}, ft(xOpt) = {:.6g}\n".format(fXOpt, ftXOpt)

        # plt.plot(N_vals, errors, label=gridToName(grid.getType(),p))
        # plt.loglog(N_vals, errors, basex=2, basey=2, label=gridToName(grid.getType(), p))
            
        # ax.plot(N_vals, errors)
        
        ax.plot(N_vals, errors, 'o', N_vals, errors, label=gridToName(grid.getType(),p))
        
    # ax.set_xscale('symlog')
    ax.set_yscale('symlog')
    # ax.xlabel('Grid Points')
    # ax.ylabel('Error')
    ax.legend(loc=3)
    plt.title("Optimization-Error Schwefel-3-1-Funktion reg-" + str(f.int_grid_level) + "-grid") 
    # plt.title("Optimization-Error Schwefel-4D-Quad-Funktion reg-3-grid") 
    plt.show()
   
def genz_error(sol):
    a = (2,2)
    exact_sol = (np.cos(a[0] + a[1]) - np.cos(a[0]) - np.cos(a[1])+1)/(a[0]*a[1])
    return abs(exact_sol-sol)

def eggcrate_error(sol):
    exact_sol = -(5/6.0)*(3*np.sin(10)-50) 
    return abs(sol - exact_sol)
    
def l2_error(d, f, f_sg):
  n = 10000
  s = 0
  for i in range(n):
    x = pysgpp.DataVector(d)
    for j in range(d):
      x[j] = random.uniform(0,1)
    s += (f_sg(x)-f(x))**2
  return (float(s)/n)**0.5

def l2_det_error(f, f_sg):
  resolution = 100
  X = np.linspace(0, 1, resolution, endpoint=True)
  s = 0
  for i in X:
    for j in X:
      vec_x = pysgpp.DataVector(2)
      vec_x[0] = i
      vec_x[1] = j
      s += (f_sg(vec_x) - f(vec_x))**2
  return (float(s)/resolution**2)**0.5
  
def integrate():
    d = 2
    grids = (pysgpp.Grid.createLinearGrid(d),
             pysgpp.Grid.createPolyGrid(d,3),
             pysgpp.Grid.createBsplineGrid(d, 3),
             # pysgpp.Grid.createLinearBoundaryGrid(d, 1),
             pysgpp.Grid.createPolyBoundaryGrid(d,3),
             pysgpp.Grid.createBsplineBoundaryGrid(d, 3),
             # pysgpp.Grid.createModLinearGrid(d),
             pysgpp.Grid.createModPolyGrid(d,3),
             pysgpp.Grid.createModBsplineGrid(d,3),
             pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
             pysgpp.Grid.createModBsplineClenshawCurtisGrid(d, 3),
             pysgpp.Grid.createModFundamentalSplineGrid(d, 3),
             pysgpp.Grid.createFundamentalSplineGrid(d, 3),
             pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 1))
    
    for grid in grids:
        #p  wird manuell gesetzt weil grid.getDegree() buggy 
        p  = 3
        if grid == grids[-1]:
            p = 1
        printLine()
        print gridToName(grid.getType(), p)
        printLine()
        try:
            hierarch = pysgpp.createOperationHierarchisation(grid)
        except:
            hierarch = pysgpp.createOperationMultipleHierarchisation(grid)
            
        # system = pysgpp.OptHierarchisationSLE(grid)
        # solver = pysgpp.OptUMFPACK()
        num_gridPoints = []
        errors = []
        try: 
            opeval = pysgpp.createOperationEval(grid)
            print("using non-naive eval")
        except:
            opeval = pysgpp.createOperationNaiveEval(grid)
            print("using naive eval")
            
        for l in range(1,9):
            grid.getStorage().clear()
            gridGen = grid.getGenerator().regular(l)
            gridStorage = grid.getStorage()
            alpha = pysgpp.DataVector(gridStorage.getSize())
            
            for i in range(gridStorage.getSize()):
                gp = gridStorage.getPoint(i)
                x = [gridStorage.getCoordinate(gp, 0), gridStorage.getCoordinate(gp ,1)]
                # alpha[i] = easom(x)
                # alpha[i] = eggcrate(x)
                alpha[i] = oscill_genz(d,x)
                # alpha[i] = target_func(x)
                # alpha[i] = michalewicz(x)
                # alpha[i] = 1
                # alpha[i] = sin_kuppel(x)
                # print(alpha[i])
            hierarch.doHierarchisation(alpha)
            # b = pysgpp.DataVector(alpha)
            # solver.solve(system, b, alpha)
  
            f_sg = lambda point: opeval.eval(alpha, point)
            # errors.append(l2_error(d, lambda point: oscill_genz(d, point), f_sg) )
            
            # errors.append(l2_det_error(lambda point: oscill_genz(d, point), f_sg) )
            # errors.append(l2_error(d, michalewicz, f_sg) )
            # errors.append(l2_error(d, easom, f_sg) )
            # errors.append(l2_error(d, eggcrate, f_sg) )
            # errors.append(l2_error(d, sin_kuppel, f_sg) )
            # try:
            opQ = pysgpp.createOperationQuadrature(grid)
            sol = opQ.doQuadrature(alpha) 
            errors.append(genz_error(sol))
            # errors.append(abs(sol - (-0.0000476373))) #easom
            # errors.append(eggcrate_error(sol)) #eggcrate
            # errors.append(abs(sol - 4/np.pi**2)) #sin_kuppel
            # except:
            # errors.append(l2_error(d, lambda point: oscill_genz(d, point), f_sg) )
            # print errors[-1]
            num_gridPoints.append(gridStorage.getSize())
            
        
        plt.loglog(num_gridPoints, errors, basex=10, basey=10, label=gridToName(grid.getType(), p))
        # plt.plot(range(1,9), errors, label=gridToName(grid.getType()))
    plt.xlabel('#Gitterpunkte')
    plt.ylabel('Fehler')
    plt.legend(loc=3)
    plt.title("Quadratur-Fehler Genz-Funktion") 
    plt.show()
    # tikz_save('mytikz.tex');

def interpolation_error():
    l = 7
    d = 2
    grids = (# pysgpp.Grid.createLinearGrid(d),
             # pysgpp.Grid.createPolyGrid(d,3),
             # pysgpp.Grid.createBsplineGrid(d, 3),
             # pysgpp.Grid.createLinearBoundaryGrid(d, 1),
             # pysgpp.Grid.createPolyBoundaryGrid(d,3),
             # pysgpp.Grid.createBsplineBoundaryGrid(d, 3),
             # pysgpp.Grid.createModLinearGrid(d),
             # pysgpp.Grid.createModPolyGrid(d,3),
             # pysgpp.Grid.createModBsplineGrid(d,3),
             # pysgpp.Grid.createLinearClenshawCurtisGrid(d),
             # pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
             # pysgpp.Grid.createModBsplineClenshawCurtisGrid(d, 3),
             # pysgpp.Grid.createModFundamentalSplineGrid(d, 3),
             pysgpp.Grid.createFundamentalSplineGrid(d, 3),
             pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 1))

    for grid in grids:
       p = 3
       if grid == grids[-1]:
           p = 1
       printLine()
       print gridToName(grid.getType(), p)
       printLine()

       grid.getGenerator().regular(l) 
       gridStorage = grid.getStorage()
       try:
           hierarch = pysgpp.createOperationHierarchisation(grid)
       except:
           hierarch = pysgpp.createOperationMultipleHierarchisation(grid)
          
       alpha = pysgpp.DataVector(gridStorage.getSize())
       nodalValues = pysgpp.DataVector(gridStorage.getSize())
       for i in range(gridStorage.getSize()):
           gp = gridStorage.getPoint(i)
           x = (gridStorage.getCoordinate(gp, 0), gridStorage.getCoordinate(gp, 1))
           alpha[i] = oscill_genz(d,x)
           # alpha[i] = sin_kuppel(x)
           # nodalValues[i] = sin_kuppel(x)
           nodalValues[i] = oscill_genz(d, x)
       hierarch.doHierarchisation(alpha)
       
       try:
           opeval = pysgpp.createOperationEval(grid)
       except: 
           opeval = pysgpp.createOperationNaiveEval(grid)
       
       f_sg = lambda point: opeval.eval(alpha, pysgpp.DataVector(point))
       
       for i in range(gridStorage.getSize()):
           gp = gridStorage.getPoint(i)
           x = [gridStorage.getCoordinate(gp, 0), gridStorage.getCoordinate(gp ,1)]
           # if(abs(sin_kuppel(x) - f_sg(x)) >= 10**-14): 
           if(abs(oscill_genz(d,x) - f_sg(x)) >= 10**-11): 
               print x
               print oscill_genz(d,x)
               # print sin_kuppel(x)
               print f_sg(x)
               print oscill_genz(d, x) - f_sg(x)
               # print sin_kuppel(x) - f_sg(x)
       resolution = 256
       X = np.linspace(0, 1, resolution, endpoint=True)
       Y = np.linspace(0, 1, resolution, endpoint=True)
       X, Y = np.meshgrid(X, Y)
       Z = np.zeros((len(X), len(X)))
       s = 0
       for i in range(len(X)):
           for j in range(len(X)):
               vec_x = pysgpp.DataVector(d)
               vec_x[0] = X[i,j] 
               vec_x[1] = Y[i,j] 
               Z[i,j] = oscill_genz(d,[X[i,j],Y[i,j]]) - f_sg(vec_x)
               # Z[i,j] = sin_kuppel(vec_x) - f_sg(vec_x)
               s += (oscill_genz(d,[X[i,j],Y[i,j]]) - f_sg(vec_x))**2
               # s += (sin_kuppel([X[i,j],Y[i,j]]) - f_sg(vec_x))**2
               # Z[i,j] = eggcrate([X[i,j],Y[i,j]]) - f_sg(vec_x)
       # opQ = pysgpp.createOperationQuadrature(grid)
       # sol = opQ.doQuadrature(alpha)
       s = (s/resolution**2)**0.5 
       fig = plt.figure()
       ax = fig.add_subplot(111, projection='3d')
       ax.plot_surface(X, Y, Z, rstride=4, cstride=4, cmap=cm.coolwarm)
       # ax.plot_wireframe(X, Y, Z)
       # plt.title("Punktweiser-Fehler Sin_Kuppel " + gridToName(grid.getType(), p) + " Level: 7 Quadratur-Fehler: "  + str(abs(sol - 4/np.pi**2)))
       plt.title("Punktweiser-Fehler Genz " + gridToName(grid.getType(), p) + " Level: 7 L2-Fehler: "  + str(s))
       plt.show()
        

def plot_genz():
  f = ExampleFunction()
  # d = f.getNumberOfParameters()
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  X = np.arange(0, 1.02, 1/50.0) 
  Y = np.arange(0, 1.02, 1/50.0)
  X, Y = np.meshgrid(X, Y)
  # Z = oscill_genz(2, [X,Y])
  # Z = schwefel([X,Y])
  # Z = sin_kuppel([X,Y])
  # Z = griewank([X,Y])
  
  Z = np.zeros((len(X),len(Y)))
  for i in range(len(X)):
      for j in range(len(X)):
          vec_x = pysgpp.DataVector(2)
          vec_x[0] = X[i,j]
          vec_x[1] = Y[i,j] 
          # Z[i,j] = f.eval(vec_x)
          # Z[i] = expeceted_revenue(vec_x[0])
          # Z[i,j] = schwefel(vec_x)
          # Z[i,j] = griewank(vec_x)
          Z[i,j] = sin_spitze(vec_x)
  
  ax.plot_surface(X, Y, Z, rstride=4, cstride=4, cmap=cm.coolwarm)
  plt.show()

def plot_1d():
  f = ExampleFunction()
  d = 1
  fig = plt.figure()
  X = np.arange(0, 1, 1/120.0) 
  Z = np.zeros(len(X))
  for i in range(len(X)):
    Z[i] = f.eval([X[i]])
    # Z[i] = -expected_revenue(X[i])
    
  plt.plot(X, Z)
  plt.show()


def int_sin(x):
  return x[0]*x[1]*np.sin(np.pi*x[0])*np.sin(np.pi.x[1])/np.pi

def test():
  grid = pysgpp.Grid.createFundamentalSplineGrid(1,3)
  grid.getGenerator().regular(5)
  gridStorage = grid.getStorage()
  print gridStorage.getSize()
  b = pysgpp.SFundamentalSplineBase(3)
  xs = np.linspace(0, 1, 200)
  for i in range(gridStorage.getSize()):
    gp = gridStorage.getPoint(i)
    x = gridStorage.getCoordinate(gp, 0)
    # print x
  i = 0
  for x in np.arange(1/64.0, 1 + 1/64.0, 1/32.0):
    print "------"
    print i
    print "x=" + str(x)
    print b.eval(5,1,x)
    i += 1
  plt.plot(xs, [b.eval(5,1, x) for x in xs])
  plt.show()

# test() 
# plot_1d()
# plot_genz()
# integrate()
optimize()
# interpolation_error()

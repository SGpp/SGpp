#!/usr/bin/python
import argparse
import pysgpp
import math
import sys
import numpy as np
import random
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from matplotlib2tikz import save as tikz_save
from pysgpp.extensions.datadriven.uq.plot import plot2d
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import checkInterpolation

class interpolation_function():
  def __init__(self, d, f):
    self.f = f
    self.d = d
    self.grid = pysgpp.Grid.createBsplineGrid(d, 3)
    self.gridStorage = self.grid.getStorage()
    # self.grid.getGenerator().regular(6)
    self.hierarch = pysgpp.createOperationMultipleHierarchisation(self.grid)
    self.opeval = pysgpp.createOperationNaiveEval(self.grid)
    self.alpha = pysgpp.DataVector(self.gridStorage.getSize())

  def create_interpolation(self, grid_lvl):
    global eval_count
    self.gridStorage.clear()
    self.grid.getGenerator().regular(grid_lvl)
    self.alpha = pysgpp.DataVector(self.gridStorage.getSize())
    for i in range(self.gridStorage.getSize()):
      gp = self.gridStorage.getPoint(i)
      x = [self.gridStorage.getCoordinate(gp, j) for j in range(self.d)]
      self.alpha[i] = self.f(x)
      eval_count += 1
    self.hierarch.doHierarchisation(self.alpha)

  def __call__(self, x):
    return self.opeval.eval(self.alpha, pysgpp.DataVector(x))


class opt_function(pysgpp.OptScalarFunction):

  def __init__(self, d):
    super(opt_function, self).__init__(d)

  def eval(self, x):
    global func, eval_count
    d = self.getNumberOfParameters()
    list_x = [x[i] for i in range(d)]
    eval_count += 1
    return func(list_x)

class opt_quad_function(pysgpp.OptScalarFunction):

  def __init__(self, d, int_d, int_grid_level):
    super(opt_quad_function, self).__init__(d)
    self.int_d = int_d
    self.int_grid_level = int_grid_level
    self.grid = pysgpp.Grid.createModBsplineClenshawCurtisGrid(self.int_d, 3)
    # self.grid = pysgpp.Grid.createBsplineGrid(self.int_d, 3)
    # self.grid = pysgpp.Grid.createLinearGrid(self.int_d)
    # self.grid = pysgpp.Grid.createModBsplineGrid(self.int_d, 3)
    self.grid.getGenerator().regular(int_grid_level)
    self.gridStorage = self.grid.getStorage()
    self.hierarch = pysgpp.createOperationMultipleHierarchisation(self.grid)
    self.quadop = pysgpp.createOperationQuadrature(self.grid)

  def eval_high_precision_quad(self, x, int_lvl):
    global func
    new_grid = pysgpp.Grid.createModBsplineClenshawCurtisGrid(int_lvl, 3)
    new_grid.getGenerator().regular(int_lvl)
    new_gridStorage = new_grid.getStorage()
    new_hierarch = pysgpp.createOperationMultipleHierarchisation(new_grid)
    new_quadop = pysgpp.createOperationQuadrature(new_grid)
    list_x = [x[i] for i in range(self.getNumberOfParameters())]
    x = list_x
    alpha = pysgpp.DataVector(new_gridStorage.getSize())
    for i in range(new_gridStorage.getSize()):
      gp = new_gridStorage.getPoint(i)
      y = [new_gridStorage.getCoordinate(gp, j) for j in range(self.int_d)]
      x_y = list(x + y)
      alpha[i] = func(x_y)
    pysgpp.OptPrinter.getInstance().setVerbosity(-2)
    new_hierarch.doHierarchisation(alpha)
    pysgpp.OptPrinter.getInstance().setVerbosity(2)
    res =  new_quadop.doQuadrature(alpha)
    # print res - schwefel(x)
    return res

  def eval(self, x):
    d = self.getNumberOfParameters()
    list_x = [x[i] for i in range(d)]
    return self.target_func(list_x)

  def target_func(self,x):
    global func, eval_count
    alpha = pysgpp.DataVector(self.gridStorage.getSize())
    for i in range(self.gridStorage.getSize()):
      gp = self.gridStorage.getPoint(i)
      y = [self.gridStorage.getCoordinate(gp, j) for j in range(self.int_d)]
      x_y = list(x + y)
      alpha[i] = func(x_y)
      eval_count += 1
    pysgpp.OptPrinter.getInstance().setVerbosity(-2)
    self.hierarch.doHierarchisation(alpha)
    pysgpp.OptPrinter.getInstance().setVerbosity(2)
    res =  self.quadop.doQuadrature(alpha)
    print res + sin_kuppel(x)
    # print res - schwefel(x)
    return res

def revenue(x):
  scaling = 100
  z = x[0]
  s = 10.0
  v = 20.0
  my = 50.0
  D = x[1]
  r = 7.0
  c = 4.0
  z *= scaling
  a = s + z
  D = (D * 2*v) + (my-v)
  assert s <= a
  return (r - c) * a + c*s - (r - 0.5*c) * max(a - D, 0)

def expected_revenue(x):
  scaling = 100
  z = x[0]
  s = 10.0
  v = 20.0
  my = 50.0
  # scaling = 2*v
  r = 7.0
  c = 4.0
  z *= scaling
  a = s + z
  assert s <= a
  expected = 0
  if (my - v < a and a <= my + v):
    expected = ((a - my + v)**2)/(4*v)
  elif(a > my + v):
    expected = a - my
  opt = (my + (v*(2*r - 3*c))/(2*r - c) - s) # / scaling
  return (r - c) * a + c*s - (r - 0.5*c) * expected

def sin_kuppel(x):
  prod = 1
  for i in range(len(x)):
    prod *= np.sin(np.pi*x[i])
  return prod

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

def ackley(x):
  x[0] -= 1
  d = len(x)
  return 20+np.e - 20*np.exp(-(np.linalg.norm(x)/(5*d**0.5))) - np.exp( (1.0/d) * sum( [ np.cos(2*np.pi*x[t]) for t in range(d)] ))

def shcb(x):
  x[0] = x[0]*10 - 5
  x[1] = x[1]*10 - 5
  return x[0]**2 * (4 - 2.1*x[0]**2 + (x[0]**4)/3) + x[0]*x[1] + 4*x[1]**2*(x[1]**2 - 1)

def peak_genz(x):
    u = [0.5, 0.5, 0.5, 0.5]
    a = [12,12,12,12]
    prod = 1
    for i in range(len(x)):
        prod *= (a[i]**-2 + (x[i] - u[i])**2)
    return 1.0/prod

def oscill_genz(x):
  global genz_a
  d = len(x)
  u = 0.5
  a = genz_a
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

def weighted_sin_cos(x):
  for i in range(len(x)):
    x[i] = x[i] * (8*(i+1.0))
  return -np.sin((np.pi/len(x)) * sum([np.cos(x[i]) for i  in range(len(x))] ))

def sin_sum(x):
  for i in range(len(x)):
    x[i] = x[i] * 5
  return -np.sin((np.pi/len(x)) * sum([x[i] for i  in range(len(x))] ))

def e_sin(x):
  for i in range(len(x)):
        x[i] = x[i] * 30
  return np.sin(np.pi*sum([x[i] for i  in range(len(x))])) * np.exp(1.0/(1+x[-1]))

def sin_nice(call_x):
  x = pysgpp.DataVector(call_x)
  n = 3
  m = 2
  d = len(x)
  assert m+n == d
  r_1 = 1
  r_2 = 4
  sin_prod = 1
  prod = 1
  for j in range(n):
    x[j] -= 0.01
    sin_prod *= np.sin(np.pi*x[j])
    prod *= x[j]

  # return schwefel(call_x[:n])+r_1*prod*np.sin(r_2*sum([x[i] for i in range(n, d - 1)])  + r_2*(2*np.pi*x[d-1]))
  return -sin_prod+r_1*np.sin(r_2*sum([x[i] for i in range(n, d - 1)])  + r_2*(2*np.pi*x[d-1]))*(np.sin(4*np.pi*sum([x[j] for j in range(n)])) + 1)
  # return sin_prod+r_1*np.sin(r_2*sum([x[i] for i in range(n, d - 1)])  + r_2*(2*np.pi*x[d-1]))*np.sin(4*np.pi*sum([x[j] for j in range(n)]))

def sin_sum_int(x):
  # return (25.0 * (np.sin((np.pi * (x[0] + x[1] + x[2]))/5.0) - 2.0 * np.sin((np.pi * (1.0 + x[0] + x[1] + x[2]))/5.0) + np.sin((np.pi * (2.0 + x[0] + x[1] + x[2]))/5.0)))/np.pi**2
  return (4/(5*np.pi))*np.cos(np.pi*1.25*(x[0] + x[1] + x[2]))- 0.8*np.cos(np.pi * 1.25*(x[0] + x[1] + x[2] +1))
  # return 4*(np.cos((5*np.pi*(x[0] + x[1]))/4) - np.cos((5*np.pi*(1 + x[0] + x[1]))/4))/(5*np.pi)

def sin_spitze(x):
  d = len(x)
  prod = 1.0
  for i in range(d):
    prod *= x[i]*np.sin(5*np.pi*x[i]+2*np.pi*x[d-i-1])
  return -prod

def griewank(x):
  #np.linalg.norm is buggy for pysgpp.DataVector
  d = len(x)
  prod = 1.0
  norm = 0.0
  for i in range(d):
    x[i] = x[i]*1200.0 - 607.0
    prod *= np.cos(x[i] / (i+1.0)**0.5 )
    norm += x[i]**2
  norm = norm**0.5
  return 1.0 + norm**2 / 4000.0 - prod

def beale(x):
  x[0] = x[0] * 10 - 5
  x[1] = x[1] * 10 - 5
  return (1.5 - x[0]*(1-x[1]))**2+(2.25 - x[0]*(1-x[1]**2))**2+(2.625-x[0]*(1-x[1]**3))**2

def mod_schwefel(call_x):
    x = pysgpp.DataVector(call_x)
    for i in range(len(x)):
        x[i] = x[i] * 1000 - 450
    return -sum([x[i]*(np.sin(abs(x[i])**0.5)) for i in range(len(x))])

def mod_schwefel_offset(call_x, int_d):
    x = pysgpp.DataVector(call_x)
    for i in range(len(x)):
        x[i] = x[i] * 1000 - 450
    return -sum([x[i]*(np.sin(abs(x[i])**0.5)) for i in range(len(x))]) + int_d*15.7759


def schwefel(call_x):
    x = pysgpp.DataVector(call_x)
    for i in range(len(x)):
        x[i] = x[i] * 1000 - 500
    return -sum([x[i]*np.sin(abs(x[i])**0.5) for i in range(len(x))]) + 30

def schwefel_prod(x):
    prod = 1.0
    for i in range(len(x)):
        x[i] = x[i] * 1000 - 500
        prod *= 10*x[i]*np.sin(abs(x[i])**0.5)
    return prod

def grid2d(start, end, num=50):
    """Create an 2D array where each row is a 2D coordinate.
    np.meshgrid is pretty annoying!
    """
    dom = np.linspace(start, end, num)
    X0, X1 = np.meshgrid(dom, dom)
    return np.column_stack([X0.flatten(), X1.flatten()])

def gridToName(i, p):
    types = {0 : "Linear(Grad "  + str(p) + ")",
             1 : "LinearStretched(Grad "  + str(p) + ")",
             2 : "LinearL0Boundary(Grad "  + str(p) + ")",
             3 : "LinearBoundary(Grad "  + str(p) + ")",
             4 : "LinearStretchedBoundary(Grad "  + str(p) + ")",
             5 : "LinearTruncatedBoundary(Grad "  + str(p) + ")",
             6 : "ModLinear(Grad "  + str(p) + ")",
             7 : "Poly (Grad " + str(p) + ")",
             8 : "PolyBoundary (Grad " + str(p) + ")",
             9 : "ModPoly(Grad "  + str(p) + ")",
             10 : "ModWavelet(Grad "  + str(p) + ")",
             11 : "ModB-Spline (Grad " + str(p) + ")",
             12 : "Prewavelet(Grad "  + str(p) + ")",
             13 : "SquareRoot(Grad "  + str(p) + ")",
             14 : "Periodic(Grad "  + str(p) + ")",
             15 : "LinearClenshawCurtis",
             16 : "B-Spline (Grad "  + str(p) + ")",
             17 : "B-SplineBoundary (Grad " + str(p) + ")",
             18 : "B-SplineClenshawCurtis (Grad " + str(p) + ")",
             19 : "Wavelet(Grad "  + str(p) + ")",
             20 : "WaveletBoundary(Grad "  + str(p) + ")",
             21 : "FundamentalSpline(Grad "  + str(p) + ")",
             22 : "ModFundamentalSpline(Grad "  + str(p) + ")",
             23 : "ModBsplineClenshawCurtis(Grad "  + str(p) + ")",
             24 : "LinearStencil(Grad "  + str(p) + ")",
             25 : "ModLinearStencil (Grad "  + str(p) + ")",
             26 : "PolyClenshawCurtisBoundary (Grad "  + str(p) + ")",
             27 : "PolyClenshawCurtis (Grad "  + str(p) + ")",
             28 : "LinearClenshawCurtis (Grad "  + str(p) + ")",
             29 : "ModPolyClenshawCurtis (Grad "  + str(p) + ")",
             30 : "ModLinearClenshawCurtis (Grad " + str(p) + ")"
             }

    return types[i]

def printLine():
    print("----------------------------------------" + \
          "----------------------------------------")


def optimize(f, sol=0, x_axis="gridpoints", use_3_grid=False):
    saved_error_lists = []
    saved_ec_lists = []
    saved_gp_lists = []
    # disable multi-threading
    global func, eval_count
    pysgpp.omp_set_num_threads(1) # increase output verbosity
    pysgpp.OptPrinter.getInstance().setVerbosity(2)
    # objective function
    # f = opt_quad_function(d, int_d, int_grid_level)
    # dimension of domain
    d = f.getNumberOfParameters()
    # B-spline degree
    p = 3
    # maximal number of grid points
    # adaptivity of grid generation
    gamma = 0.85
    # N_vals = [50, 100, 200, 500, 1000, 2000, 4000]
    # N_vals = [1000, 2000]
    # N_vals = [1000]
    N_vals = [1, 10, 50, 100, 200, 500, 1000, 2000]#, 500, 1000] # 2000]
    # N_vals = [5, 15, 20, 30, 40]
    loop_list = []
    if(use_3_grid):
      interpolf = interpolation_function(d + f.int_d, func)
      print("using 3-grid method")
      loop_list= range(1,8)
    else:
      loop_list = N_vals
      grids = [
             # (pysgpp.Grid.createBsplineGrid(d, 3), 3),
             # (pysgpp.Grid.createBsplineBoundaryGrid(d, 3), 3),
             (pysgpp.Grid.createModBsplineGrid(d, 3), 3),
             # (pysgpp.Grid.createModBsplineClenshawCurtisGrid(d, 3),3),
             # (pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3), 3)
    ]
    print("starting optimization with {} different grids, plotting based on {}, up to {} Ritter-Novak points".format(len(grids), x_axis, N_vals[-1]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for grid, p in grids:
        gridpoint_list = []
        ec_list = []
        printLine()
        printLine()
        errors = []
        for loop_item in loop_list:
            eval_count = 0
            if(use_3_grid):
              grid_lvl = loop_item
              interpolf.create_interpolation(grid_lvl)
              func = interpolf
              N = 3000
            else:
              N = loop_item
            grid.getStorage().clear()
            gridGen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma)

            # #############################################################################
            # GRID GENERATION
            # #############################################################################

            printLine()
            print("Generating {}".format(gridToName(grid.getType(),p)))


            if not gridGen.generate():
                print "Grid generation failed, exiting."
                sys.exit(1)
            gridpoint_list.append(grid.getStorage().getSize())
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
            ftHessian = pysgpp.OptInterpolantScalarFunctionHessian(grid, coeffs)

            gradientMethod = pysgpp.OptAdaptiveGradientDescent(ft, ftGradient)
            newtonsMethod = pysgpp.OptAdaptiveNewton(ft, ftHessian)
            # differentialEvolution = pysgpp.OptDifferentialEvolution(ft)

            optimizers = [pysgpp.OptMultiStart(gradientMethod, 100000, 300),
                          pysgpp.OptMultiStart(newtonsMethod, 100000, 300)]
                          # differentialEvolution]

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

            xOpt, fXOpt, ftXOpt = x0, fX0, ftX0

            for optimizer in optimizers:
              optimizer.optimize()
              xOptCur = optimizer.getOptimalPoint()
              fXOptCur = f.eval(xOptCur)
              ftXOptCur = optimizer.getOptimalValue()

              if fXOptCur < fXOpt:
                xOpt, fXOpt, ftXOpt = xOptCur, fXOptCur, ftXOptCur

            # errors.append(1 + fXOpt) # 80.00256374
            ec_list.append(eval_count)
            # use real funct for error calculation only really makes sense for schwefel function
            # errors.append(schwefel(xOpt)  - sol)
            # errors.append(mod_schwefel_offset(xOpt, int_d)  - sol)
            if isinstance(f, opt_quad_function):
              # pass
              # errors.append(f.eval_high_precision_quad(xOpt, 5) - sol)
              x_temp = [xOpt[i] - 0.01 for i in range(d)]
              errors.append(1-sin_kuppel(x_temp))
            elif isinstance(f, opt_function):
              errors.append(fXOpt - sol)
            # errors.append(abs(fXOpt - sol))
            print "\nxOpt = {}".format(xOpt)
            print "f(xOpt) = {:.13f}, ft(xOpt) = {:.13f}".format(fXOpt, ftXOpt)
            print "sol = {:.8f}, error = {:.2e}\n".format(sol, errors[-1])
            # print "mod_schwefel_offset(xOpt)= {}".format(mod_schwefel_offset(xOpt, int_d))
            print "Eval Count: {}".format(eval_count)

        color = "k"
        degree_based = False
        if "ModBsplineClenshaw" in gridToName(grid.getType(), p):
            color = "g"
        elif "Boundary" in gridToName(grid.getType(), p):
            color = "r"
        elif "Clenshaw" in gridToName(grid.getType(), p):
            color = "k"
        elif "Mod" in gridToName(grid.getType(), p):
            color = "c"
        elif "B-Spline" in gridToName(grid.getType(), p):
            color = "b"
        linestyle = color + ".-"
        saved_ec_lists.append(list(ec_list))
        saved_error_lists.append(list(errors))
        saved_gp_lists.append(list(gridpoint_list))
        if(x_axis == "gridpoints"):
          ax.plot(gridpoint_list, errors, linestyle, label=gridToName(grid.getType(),p))
          plt.xlabel('\#Gitterpunkte')
        elif(x_axis == "evaluations"):
          ax.plot(ec_list, errors, linestyle, label=gridToName(grid.getType(),p))
          plt.xlabel('Evaluationen')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.ylabel(r'$v(\vec{x}^{*}) - v(\vec{x}_{opt})$')
    print "errors:"
    print saved_error_lists
    print "gridpoint list:"
    print saved_gp_lists
    print "ec list:"
    print saved_ec_lists

    # plt.ylabel('$v(\vec{x}^{*}) - v(\vec{x}_{\text{opt}})$')
    # ax.legend(loc=0)

    if isinstance(f, opt_function):
        # plt.title("Optimization-Errorotation=0, fontsize=16)
        save_name = 'opt-{}-{}'.format(func_name, d)
    elif isinstance(f, opt_quad_function):
        # plt.title("Optimization-Error {}-{}-{}-Funktion reg-{}-grid".format(func_name, d, f.int_d, f.int_grid_level))
        save_name= 'opt-quad-{}-{}-{}-{}'.format(func_name, d, f.int_d, f.int_grid_level)
    tikz_save(save_name + ".tex", figureheight = '\\figureheight', figurewidth = '\\figurewidth')
    plt.show()
    # plt.savefig(save_name + ".png")

def genz_error():
    a = genz_a
    d = len(a)
    if d == 2:
      exact_sol = (np.cos(a[0] + a[1]) - np.cos(a[0]) - np.cos(a[1])+1)/(a[0]*a[1])
    elif all([i == 2 for i in a]):
      return np.sin(1)**d * (-np.cos(d))
    else:
      print "no exact solution for this genz function"
    return exact_sol

def eggcrate_error():
    exact_sol = -(5/6.0)*(3*np.sin(10)-50)
    return exact_sol

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

def integrate(d, error_type = "l2", quad_error_sol = 0):
    mixedGrids = [ #  (pysgpp.Grid.createLinearGrid(d), 1),
             # (pysgpp.Grid.createPolyGrid(d,3), 3),
             # (pysgpp.Grid.createBsplineGrid(d, 3), 3) ,
             # (pysgpp.Grid.createBsplineGrid(d, 5), 5) ,
             # (pysgpp.Grid.createBsplineGrid(d, 7), 7),
             # (pysgpp.Grid.createLinearBoundaryGrid(d, 1), 1),
             (pysgpp.Grid.createPolyClenshawCurtisGrid(d, 3), 3),
             (pysgpp.Grid.createModPolyClenshawCurtisGrid(d, 3), 3),
             (pysgpp.Grid.createPolyBoundaryGrid(d,3), 3),
             (pysgpp.Grid.createBsplineBoundaryGrid(d, 3), 3),
             (pysgpp.Grid.createLinearClenshawCurtisGrid(d), 1),
             (pysgpp.Grid.createModLinearClenshawCurtisGrid(d), 1),
             # (pysgpp.Grid.createModLinearGrid(d), 1),
             # (pysgpp.Grid.createModPolyGrid(d,3), 3),
             (pysgpp.Grid.createModBsplineGrid(d,3), 3),
             (pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3), 3),
             (pysgpp.Grid.createModBsplineClenshawCurtisGrid(d, 3), 3),
             (pysgpp.Grid.createModFundamentalSplineGrid(d, 3), 3),
             # (pysgpp.Grid.createFundamentalSplineGrid(d, 3), 3 ),
             # (pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 1), 1),
    ]
    plot_mode = False# True -> bsplineGrids Flase -> mixedGrids
    # bsplineGrids =[ # (pysgpp.Grid.createBsplineGrid(d, 3), 3),
    #                # (pysgpp.Grid.createBsplineGrid(d, 5), 5),
    #                # (pysgpp.Grid.createBsplineGrid(d, 7), 7),
    #                (pysgpp.Grid.createBsplineBoundaryGrid(d, 3), 3),
    #                (pysgpp.Grid.createBsplineBoundaryGrid(d, 5), 5),
    #                (pysgpp.Grid.createBsplineBoundaryGrid(d, 7), 7),
    #                (pysgpp.Grid.createModBsplineGrid(d, 3), 3),
    #                (pysgpp.Grid.createModBsplineGrid(d, 5), 5),
    #                (pysgpp.Grid.createModBsplineGrid(d, 7), 7),
    #                (pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3), 3),
    #                (pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 5), 5),
    #                (pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 7), 7),

    # ]


    for grid, p in mixedGrids:
    # for grid, p in bsplineGrids:
        #p  wird manuell gesetzt weil grid.getDegree() buggy
        printLine()
        print gridToName(grid.getType(), p)
        printLine()
        try:
            hierarch = pysgpp.createOperationHierarchisation(grid)
        except:
            hierarch = pysgpp.createOperationMultipleHierarchisation(grid)

        num_gridPoints = []
        errors = []
        if error_type == "l2":
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
                x = list([gridStorage.getCoordinate(gp, c) for c in range(d)])
                alpha[i] = func(x)
            hierarch.doHierarchisation(alpha)

            f_sg = lambda point: opeval.eval(alpha, point)

            if error_type == "quad":
                opQ = pysgpp.createOperationQuadrature(grid)
                sol = opQ.doQuadrature(alpha)
                errors.append(abs(quad_error_sol - sol))
            elif error_type == "l2":
                errors.append(l2_error(d, func, f_sg))
            num_gridPoints.append(gridStorage.getSize())

        linestyle = ".-"
        color = "k"
        degree_based = False
        if not plot_mode:
            if "Mod" in gridToName(grid.getType(), p):
                linestyle = ".--"
            elif "Boundary" in gridToName(grid.getType(), p):
                linestyle = ".-"

            if "Clenshaw" in gridToName(grid.getType(), p):
                color = "k"
            elif "B-Spline" in gridToName(grid.getType(), p):
                color = "b"
            elif "Poly" in gridToName(grid.getType(), p):
                color = "g"
            elif "Linear" in gridToName(grid.getType(), p):
                color = "r"
            elif "Fundamental" in gridToName(grid.getType(), p):
                color = "c"
        else:
            if p == 3:
                color = "k"
            elif p == 5:
                color = "b"
            elif p == 7:
                color = "r"

            if "Mod" in gridToName(grid.getType(), p):
                linestyle = ".--"
            elif "Boundary" in gridToName(grid.getType(), p):
                linestyle = ".:"

        linestyle = color + linestyle
        plt.loglog(num_gridPoints, errors, linestyle, basex=10, basey=10, label=gridToName(grid.getType(), p) )
    plt.xlabel(r'\#Gitterpunkte')
    # plt.legend(loc=3)
    if(error_type == "quad"):
        save_title = "quad_{}.tex".format(func_name)
        plt.ylabel(r'$|\epsilon|$', fontsize=16)
        # plt.title("Quadratur-Fehler {}-Funktion".format(func_name))
    elif(error_type == "l2"):
        # plt.title("L2-Fehler {}-Funktion".format(func_name))
        plt.ylabel(r'$L_{2}$-Fehler', fontsize=16)
        save_title = "l2_{}.tex".format(func_name)

    tikz_save(save_title, figureheight = '\\figureheight', figurewidth = '\\figurewidth');
    plt.show()

def interpolation_error():
    l = 7
    d = 2
    grids = [# pysgpp.Grid.createLinearGrid(d),
             # pysgpp.Grid.createPolyGrid(d,3),
             # pysgpp.Grid.createBsplineGrid(d, 3),
             # pysgpp.Grid.createLinearBoundaryGrid(d, 1),
             pysgpp.Grid.createPolyBoundaryGrid(d,3),
             pysgpp.Grid.createBsplineBoundaryGrid(d, 3),
             # pysgpp.Grid.createModLinearGrid(d),
             # pysgpp.Grid.createModPolyGrid(d,3),
             pysgpp.Grid.createModBsplineGrid(d,3),
             # pysgpp.Grid.createLinearClenshawCurtisGrid(d),
             pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3),
             # pysgpp.Grid.createModBsplineClenshawCurtisGrid(d, 3),
             # pysgpp.Grid.createModFundamentalSplineGrid(d, 3),
             # pysgpp.Grid.createFundamentalSplineGrid(d, 3),
             # pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 1))
        ]
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
       for i in range(gridStorage.getSize()):
           gp = gridStorage.getPoint(i)
           x = [gridStorage.getCoordinate(gp, 0), gridStorage.getCoordinate(gp, 1)]
           alpha[i] = func(x)
       hierarch.doHierarchisation(alpha)

       try:
           opeval = pysgpp.createOperationEval(grid)
       except:
           opeval = pysgpp.createOperationNaiveEval(grid)

       f_sg = lambda point: opeval.eval(alpha, pysgpp.DataVector(point))

       resolution = 150
       for x in np.linspace(0, 1, 200):
           print f_sg([x,0]) - func([x,0])
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
               Z[i,j] = func(vec_x) - f_sg(vec_x)
               s += (func([X[i,j],Y[i,j]]) - f_sg(vec_x))**2
       # opQ = pysgpp.createOperationQuadrature(grid)
       # sol = opQ.doQuadrature(alpha)
       s = (s/resolution**2)**0.5
       fig = plt.figure()
       ax = fig.add_subplot(111, projection='3d')
       ax.plot_surface(X, Y, Z, rstride=4, cstride=4, cmap=cm.coolwarm)
       z_formatter = ticker.ScalarFormatter()
       z_formatter.set_powerlimits((-1,1))
       ax.zaxis.set_major_formatter(z_formatter)
       ax.set_xlabel(r'x$_{1}$', fontsize=16)
       ax.set_ylabel(r'x$_{2}$', fontsize=16)
       # ax.zaxis.set_rotate_label(False)
       # ax.set_zlabel(r'Fehler', rotation=0, fontsize=16)
       ax.set_zlabel(r'$q(\vec{x}) - \widetilde{q}(\vec{x})$', fontsize=16)
       # ax.plot_wireframe(X, Y, Z)
       # plt.title("Punktweiser-Fehler Sin_Kuppel " + gridToName(grid.getType(), p) + " Level: 7 Quadratur-Fehler: "  + str(abs(sol - 4/np.pi**2)))
       # plt.title("Punktweiser-Fehler {} {} Level: 7 L2-Fehler: {}".format(func_name, gridToName(grid.getType(), p), s))
       # tikz_save("pointwise_error_{}_{}.tex".format(func_name, gridToName(grid.getType(), p) ), figureheight = '\\figureheight', figurewidth = '\\figurewidth');
       plt.show()


def plot_func(f):
  # f = ExampleFunction()
  # d = f.getNumberOfParameters()
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  X = np.linspace(0, 1, 50, endpoint=True)
  Y = np.linspace(0, 1, 50, endpoint=True)
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

          Z[i,j] = f(vec_x)
          # Z[i] = expeceted_revenue(vec_x[0])
          # Z[i,j] = schwefel(vec_x)
          # Z[i,j] = griewank(vec_x)
          # Z[i,j] = sin_spitze(vec_x)
          # Z[i,j] = func(vec_x)
          # Z[i,j] = func([X[i,j], Y[i,j]])

          # Z[i,j] = f.eval([X[i,j], Y[i,j]])
  ax.set_xlabel(r'x$_{1}$', fontsize=16)
  ax.set_ylabel(r'x$_{2}$', fontsize=16)
  # ax.zaxis.set_rotate_label(False)
  ax.set_zlabel(r'$q(\vec{x})$', fontsize=16)

  ax.plot_surface(X, Y, Z, rstride=4, cstride=4, cmap=cm.coolwarm)
  # tikz_save("{}-funktion.tex".format(func_name), figureheight = '\\figureheight', figurewidth = '\\figurewidth');
  plt.show()

def plot_1d():
  # f = opt_quad_function()
  # d = 1
  # fig = plt.figure()
  # X = np.arange(0, 1, 1/120.0)
  X = np.arange(0, 1, 1/100.0)
  # for f in functions:
  Z = np.zeros(len(X))
  for i in range(len(X)):
    # Z[i] = -func([X[i]])
    Z[i] = func(X[i])
    # Z[i] = -expected_revenue(X[i])
  plt.xlabel(r"z")
  plt.ylabel(r"v(z)", rotation=0)
  plt.plot(X, Z)
  plt.show()

func_dic = {"genz" : oscill_genz,
            "ackley" : ackley,
            "sin_kuppel" : sin_kuppel,
            "sin_sum" : sin_sum,
            "sin_sum_int" : sin_sum_int,
            "sin_nice" : sin_nice,
            "revenue" : revenue,
            "expected_revenue" : lambda point: -expected_revenue(point),
            "eggholder" : eggholder,
            "eggcrate" : eggcrate,
            "schwefel" : schwefel,
            "mod_schwefel" : mod_schwefel,
            "schwefel_prod" : schwefel_prod,
            "shcb" : shcb,
            "easom" : easom,
            "griewank" : griewank,
            "beale" : beale,
            "michalewicz" : michalewicz,
            "sin_cos" : weighted_sin_cos,
            "peak_genz" : peak_genz,
            "e_sin" : e_sin}

def quad_sol(function_name, d):
  if function_name == "sin_kuppel":
    return (2/np.pi)**d
  elif function_name == "eggcrate":
    return eggcrate_error()
  elif function_name == "mod_schwefel":
    return 15.7759*d
  elif function_name == "schwefel":
    return 30.0
  elif function_name == "peak_genz":
    return np.pi**2
  elif function_name == "easom":
    return -0.0000476373
  elif function_name == "genz":
    return genz_error(d)

def opt_sol(function_name, d):
  if(function_name == "schwefel"):
    return d*(-418.982887272433743) + 30
  if(function_name == "mod_schwefel"):
    return d*(-418.982887272433743)+ 2*15.7759
    # return mod_schwefel([0.924928405301430500707 for i in range(d)])+15.7759*2
  elif(function_name == "sin_nice"):
    return -0.0042070893493 #lvl 5
    return -0.6714977322570 #lvl 5
    return -0.6714977346292 #lvl 7
    # return -0.679251695848
    # return sin_nice([6.42854853427698991197e-01, 3.92701901570129685570e-01])
    # return -0.9003163161571061
  elif(function_name == "sin_nice_int"):
    return -0.9003163161571061
  elif(function_name == "expected_revenue"):
    return -expected_revenue([0.44])

# elif(args.function == "sin_sum"):
      # sol = -0.967531
  else:
    return 0

if __name__ == "__main__":
  global genz_a
  parser = argparse.ArgumentParser()
  operation_group = parser.add_mutually_exclusive_group()
  operation_group.add_argument("-oq", "--optimize_quad", help="perform optimization of a quadrature function", nargs = 4)
  operation_group.add_argument("-op", "--optimize", help="perform optimization", nargs=2)
  operation_group.add_argument("-c", "--convergence", help="create convergence plot", nargs=2)
  operation_group.add_argument("-p", "--plot", help="plot a function", nargs=1)
  operation_group.add_argument("-e", "--pointwise_error", help="plot the pointwise error of a function", action="store_true")
  parser.add_argument("function", help="function")
  parser.add_argument("-g_3", "--use_3_grid", help="create overall interpolation before optimizing, onmly makles sense with -oq", action="store_true", default=False)
  args = parser.parse_args()
  # func = schwefel_prod
  # f = opt_quad_function(3,1,1)
  # print f.eval([0.1, 0.1, 0.1] )
  # print args
  func = func_dic[args.function]
  func_name = args.function[0].upper() + args.function[1:]
  eval_count = 0
  # print "{:.16f}".format(sin_nice_int([0.5, 0.5+1/7.0]))
  if args.optimize_quad != None:
    d = int(args.optimize_quad[0])
    int_d = int(args.optimize_quad[1])
    int_level = int(args.optimize_quad[2])
    x_axis = args.optimize_quad[3] #evaluations or gridpoints
    if(x_axis != "evaluations" and x_axis != "gridpoints"):
      print "Invalid argument for x_axis"
      exit()
    f = opt_quad_function(d, int_d, int_level)
    sol = opt_sol(args.function, d)
    # plot_func(f.eval)
    # print("{:.13f}").format(f.eval([6.42854853427698991197e-01, 3.92701901570129685570e-01]))
    sol = 0
    optimize(f, sol, x_axis, args.use_3_grid)
  elif args.convergence != None:
    d = int(args.convergence[0])
    genz_a = [2 for i in range(d)]
    error_type = args.convergence[1]
    if args.convergence[1] == "quad":
      quad_error_sol = quad_sol(args.function, d)
      if quad_error_sol == None:
        quad_error_sol = 0
      integrate(d, error_type, quad_error_sol)
    else:
      integrate(d, error_type)
  elif args.plot != None:
    d = int(args.plot[0])
    genz_a = [2 for i in range(d)]
    if d == 1:
        plot_1d()
    elif d == 2:
        plot_func(func)
  elif args.pointwise_error:
    interpolation_error()
  elif args.optimize != None:
    d = int(args.optimize[0])
    x_axis = args.optimize[1] #evaluations or gridpoints
    sol = opt_sol(args.function, d)
    f = opt_function(d)
    optimize(f, sol, x_axis, args.use_3_grid)
    # plot_1d()

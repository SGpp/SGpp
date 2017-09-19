import pysgpp

import matplotlib.pyplot as plt
import numpy as np
import random
from pprint import pprint as pp
from itertools import chain, combinations
import scipy.misc as spmisc
from mpl_toolkits.mplot3d import Axes3D


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def operationwrapper(combi_operation, level):
    def operation(x):
        return combi_operation.evaluate(level, x)

    return operation


# for functions from sgpp::optimization
def getfuncwrapper(func):
    def function(x):
        return func.evalUndisplaced(x)

    return function


# get gradfunction in shape f(x)
def getgradkfunc(func, k):
    def gradk(x):
        return grad_xk(x, k, func)

    return gradk


def estimatel2Error(n, dim, gridOpEval, targetFunc):
    sum = 0
    for _ in range(n):
        point = pysgpp.DataVector(dim)
        for d in range(dim):
            point[d] = random.random()
        # print(gridOpEval(point)-targetFunc.evalUndisplaced(point))
        sum += (gridOpEval(point) - targetFunc(point)) ** 2
    sum = sum / n
    sum = np.sqrt(sum)
    return sum


## todo wrong because eval(0) is by definition 0?
def calc_error_on_gridpts(gridpoints, gridOpEval, targetFunc):
    error = 0

    for i in range(len(gridpoints)):
        x = gridpoints[i]
        error = max(gridOpEval(x) - targetFunc(x), error)
        # print(str(x[0]) + "   " + str(x[1]), abs(gridOpEval(x) - targetFunc(x)), gridOpEval(x),
        # targetFunc(x))

    return error


def calc_error_gradient(gridpoints, gridOpEval, targetFunc, dim):
    error = []

    # iterate over dimensions
    for d in range(dim):
        grad_i = getgradkfunc(targetFunc, d)
        opEvalgrad_i = getgradkfunc(gridOpEval, d)
        error.append(calc_error_on_gridpts(gridpoints, opEvalgrad_i, grad_i))
        # print("-------------------------------------------------------------")

    return error


def calc_error_mixed_gradient(gridpoints, gridOpEval, targetFunc, dim):
    mixgrad = targetFunc
    opEval_mixgrad = gridOpEval
    for d in range(dim):
        mixgrad = getgradkfunc(mixgrad, d)
        opEval_mixgrad = getgradkfunc(opEval_mixgrad, d)
    # print("-------------------------------------------------------------")
    error = calc_error_on_gridpts(gridpoints, opEval_mixgrad, mixgrad)
    return error


def displace(x, i, h):
    x_temp = pysgpp.DataVector(x)
    x_temp[i] = x[i] + h
    return x_temp


# get the derivatives for point x
def grad_x(x, function):
    h = 1e-09
    gy = np.zeros(x.shape, dtype=np.float64)

    for k in range(x.shape[0]):
        gy[k] = (function(displace(x, k, h)) -
                 function(displace(x, k, -h))) / (2 * h)

    return gy


"""
x evaluationpoint

returns the k-th gradient with finite differences method
"""


def grad_xk(x, k, function):
    h = 1e-06

    if x[k] == 0:
        gy = (function(displace(x, k, h)) - function(x)) / (h)

    elif x[k] == 1:

        gy = (function(x) - function(displace(x, k, -h))) / (h)
    else:

        gy = (function(displace(x, k, h)) - function(displace(x, k, -h))) / (2 * h)
    return gy


def calc_tangent(x, k, operation):
    grad_k = getgradkfunc(operation, k)

    m = grad_k(x)
    b = operation(x)
    return m, b


# returns a number of points to draw the tangent
def tangent_samples(x, dim, k, m, b):
    samples = 2
    x_samples = []
    for d in range(dim):
        x_samples.append(np.full(samples, x[d]))

    x_samples[k] = (np.linspace(x[k] - 0.05, x[k] + 0.05, samples))

    y = [(x_samples[k][i] - x[k]) * m + b for i in range(samples)]

    return x_samples[0], x_samples[1], y


def calc_gradient_tangent_gridpoint(gridpoints, operation, dim):
    x0_all = []
    x1_all = []
    y_all = []
    for p in range(len(gridpoints)):

        for d in range(dim):
            m, b = calc_tangent(gridpoints[p], d, operation)

            x0, x1, y = tangent_samples(gridpoints[p], 2, d, m, b)
            x0_all.append(x0)
            x1_all.append(x1)
            y_all.append(y)

    return x0_all, x1_all, y_all


# parabola between [0,1]
def f1D(x):
    return -4 * ((x[0] - 0.5) ** 2) + 1


def f1D_grad(x):
    return 4 - 8 * x[0]


def f2D(x):
    return 4 * x[0] * (1 - x[0]) * 4 * x[1] * (1 - x[1])


def f2D_test(x):
    return x[0] ** 2 + x[1] ** 2


# combi_combigrids

class CombiCombigrid1d_hermite:
    def __init__(self, function, grad_function):
        self.func = pysgpp.multiFunc(function)
        self.grad_func = pysgpp.multiFunc(grad_function)
        self.d = 1
        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)
        self.operation_zeta = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
            self.d, self.grad_func)

    def evaluate(self, level, x):
        return self.operation_psi.evaluate(level,
                                           x) + 1 * self.operation_zeta.evaluate(
            level, x)

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()



class CombiCombigrid2dHermite_without_mixed:
    def __init__(self, function):
        self.func = pysgpp.multiFunc(function)
        self.grad_0 = pysgpp.multiFunc(getgradkfunc(self.func, 0))
        self.grad_1 = pysgpp.multiFunc(getgradkfunc(self.func, 1))
        # derivative 0 then 1
        self.mixed_grad = pysgpp.multiFunc(getgradkfunc(self.grad_0, 1))

        self.d = 2
        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)
        self.operation_psi_zeta_0 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, 0, self.grad_0)
        self.operation_psi_zeta_1 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, 1, self.grad_1)
        self.operation_zeta = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
            self.d, self.mixed_grad)

    def evaluate(self, level, x):
        sum = 0
        sum += self.operation_psi.evaluate(level, x)
        sum += self.operation_psi_zeta_0.evaluate(level, x)
        sum += self.operation_psi_zeta_1.evaluate(level, x)
        # sum += self.operation_zeta.evaluate(level, x)
        return sum

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()


class CombiCombigridHermite:
    def __init__(self, function, dim):
        self.d = dim
        self.func = pysgpp.multiFunc(function)
        self.grad = []
        ##derivatives
        for i in range(self.d):
            self.grad.append(pysgpp.multiFunc(getgradkfunc(self.func, i)))

        self.mixed_grad = []
        for i in range(self.d):
            self.mixed_grad = getgradkfunc(self.mixed_grad, i)
        self.mixed_grad = pysgpp.multiFunc(self.mixed_grad)
        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)

        self.operation_zeta_psi = []
        for i in range(self.d):
            self.operation_zeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(self.d,
                                                                                          i,
                                                                                          self.grad[
                                                                                              i]))
        self.mixed_directions = list(powerset([i for i in range(self.d)]))

        self.mixed_grad = []

        self.operationzeta_psi = []
        for i in range(1, len(self.mixed_directions)):

            mixed_grad_temp = self.func
            for j in self.mixed_directions[i]:
                mixed_grad_temp = getgradkfunc(mixed_grad_temp, j)
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))

        for i in range(len(self.mixed_grad)):
            print(self.mixed_directions[i])
            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    self.mixed_directions[
                                                                                        i + 1],
                                                                                    self.mixed_grad[
                                                                                        i]))

    def evaluate(self, level, x):
        sum = 0

        sum += self.operation_psi.evaluate(level, x)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += op.evaluate(level, x)
        return sum

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()


class HierachGridBSpline:
    ##create a bsplineGrid with calculated hierarchized coefficients

    def __init__(self, dim, degree, func):
        self.dim = dim
        self.degree = degree
        self.func = function_eval_undisplaced(func)
        self.level = -1

    def bsplineGrid(self, level):

        grid = pysgpp.Grid.createBsplineGrid(self.dim, self.degree)
        self.gridStorage = grid.getStorage()
        print "dimensionality:         {}".format(self.gridStorage.getDimension())

        grid.getGenerator().regular(level)

        functionValues = pysgpp.DataVector(self.gridStorage.getSize())
        functionValues.setAll(0.0)
        for i in xrange(self.gridStorage.getSize()):
            gp = self.gridStorage.getPoint(i)
            coordinates = pysgpp.DataVector(self.dim)
            gp.getStandardCoordinates(coordinates)
            # print( coordinates[0],coordinates[1])
            functionValues[i] = self.func.evalUndisplaced(coordinates)

        print ("Hierarchizing...\n")
        coeffs = pysgpp.DataVector(len(functionValues))
        hierSLE = pysgpp.OptHierarchisationSLE(grid)
        sleSolver = pysgpp.OptAutoSLESolver()

        if not sleSolver.solve(hierSLE, functionValues, coeffs):
            print "Solving failed, exiting."

        print(self.level)

        return grid, coeffs

    def gridOperationEval(self, opEval, alpha):
        def eval(x):
            return opEval.eval(alpha, x)

        return eval

    def evaluate(self, level, x):
        if self.level != level + 1:
            self.level = level + 1
            self.grid, self.coeffs = self.bsplineGrid(self.level)
            self.opEvalBspline = pysgpp.createOperationEvalNaive(self.grid)

        return self.opEvalBspline.eval(self.coeffs, x)

    def getLevelManager(self):
        return self.LevelManager(self.gridStorage)

    class LevelManager:
        def __init__(self, gridstorage):
            self.gridStorage = gridstorage

        def numGridPoints(self):
            return self.gridStorage.getSize()


class function_eval_undisplaced:
    def __init__(self, function):
        self.function = function

    def evalUndisplaced(self, x):
        return self.function(x)


class CombiCombigrid2dLinear:
    def __init__(self, function):
        self.func = pysgpp.multiFunc(function)
        self.grad_0 = pysgpp.multiFunc(getgradkfunc(self.func, 0))
        self.grad_1 = pysgpp.multiFunc(getgradkfunc(self.func, 1))
        self.d = 2
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)
        self.operation_psi_0 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(
            self.d, 0, self.func)
        self.operation_psi_1 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(
            self.d, 1, self.func)
        self.operation_zeta_0 = pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
            self.d, 0, self.grad_0)
        self.operation_zeta_1 = pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
            self.d, 1, self.grad_1)

    def evaluate(self, level, x):
        sum = 0
        sum += self.operation_psi_0.evaluate(level, x) + \
               self.operation_zeta_0.evaluate(level, x)
        sum += self.operation_psi_1.evaluate(level, x) + \
               self.operation_zeta_1.evaluate(level, x)
        sum -= self.operation_linear.evaluate(level, x)
        return sum

    def getLevelManager(self):
        return self.operation_linear.getLevelManager()



class CombiCombigridLinear:
    def __init__(self, function, dim):
        self.func = pysgpp.multiFunc(function)
        self.grad = []
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(getgradkfunc(self.func, i)))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_psi = []
        self.operation_zeta = []
        for i in range(dim):
            self.operation_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(self.d,
                                                                                         i,
                                                                                         self.func))

        for i in range(dim):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
                    self.d, i, self.grad[i]))

    def evaluate(self, level, x):
        sum = 0

        for operation in self.operation_psi:
            sum += operation.evaluate(level, x)
        for operation in self.operation_zeta:
            sum += operation.evaluate(level, x)
        sum -= self.operation_linear.evaluate(level, x) * (self.d - 1)
        return sum

    def getLevelManager(self):
        return self.operation_linear.getLevelManager()

################FullGrids########################################################
class CombiFullGridHermite:
    def __init__(self, function, dim):
        self.d = dim
        self.func = pysgpp.multiFunc(function)
        self.grad = []
        self.level=0
        ##derivatives
        for i in range(self.d):
            self.grad.append(pysgpp.multiFunc(getgradkfunc(self.func, i)))

        self.mixed_grad = []
        for i in range(self.d):
            self.mixed_grad = getgradkfunc(self.mixed_grad, i)
        self.mixed_grad = pysgpp.multiFunc(self.mixed_grad)
        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)
        self.operation_psi = self.operation_psi.getFullGridEvaluator()

        self.operation_zeta_psi = []
        for i in range(self.d):
            self.operation_zeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                    self.d,
                    i,
                    self.grad[
                        i]))
        self.mixed_directions = list(powerset([i for i in range(self.d)]))

        self.mixed_grad = []

        self.operationzeta_psi = []
        for i in range(1, len(self.mixed_directions)):

            mixed_grad_temp = self.func
            for j in self.mixed_directions[i]:
                mixed_grad_temp = getgradkfunc(mixed_grad_temp, j)
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))

        for i in range(len(self.mixed_grad)):

            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    self.mixed_directions[
                                                                                        i + 1],
                                                                                    self.mixed_grad[
                                                                                        i]))
        #lazy full conversion to full grid operations

        for i in range(len(self.operationzeta_psi)):
            self.operationzeta_psi[i]=self.operationzeta_psi[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level=level

        sum += self.eval_full(level,x,self.operation_psi)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += self.eval_full(level,x,op)
        return sum


    def eval_full(self, level, x, fullgrid_operation):
        scalars = []
        for i in range(len(x)):

            scalars.append(pysgpp.FloatScalarVector(x[i]))

        scalars=pysgpp.FloatScalarVectorVector(scalars)
        fullgrid_operation.setParameters(scalars)

        return fullgrid_operation.eval([level for i in range(self.d)]).getValue()

    def getLevelManager(self):
        return self.LevelManager(self.operation_psi,self.level,self.d)

    class LevelManager:
        def __init__(self,operation,level,d):
            self.operation=operation
            self.level=level
            self.d=d


        #returns points with fullgrid
        def numGridPoints(self):
            return self.operation.numPoints([self.level for i in range(self.d)])
        def getAllGridPoints(self):
            return  self.operation.getGridPoints([self.level for i in range(self.d)])


class CombiFullGridLinear:
    def __init__(self, function, dim):
        self.func = pysgpp.multiFunc(function)
        self.grad = []
        self.level=0
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(getgradkfunc(self.func, i)))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_linear=self.operation_linear.getFullGridEvaluator()



        self.operation_psi = []
        self.operation_zeta = []
        for i in range(dim):
            self.operation_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(self.d,
                                                                                         i,
                                                                                         self.func))
        # lazy redefine to fullgrid
        for i in range(dim):
            self.operation_psi[i] = self.operation_psi[i].getFullGridEvaluator()

        for i in range(dim):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
                    self.d, i, self.grad[i]))

        # lazy redefine to fullgrid
        for i in range(dim):
            self.operation_zeta[i]=self.operation_zeta[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level=level

        for operation in self.operation_psi:
            sum += self.eval_full(level,x,operation)
        for operation in self.operation_zeta:
            sum += self.eval_full(level,x,operation)
        sum -= self.eval_full(level,x,self.operation_linear) * (self.d - 1)
        return sum

    def eval_full(self, level, x, fullgrid_operation):
        scalars = []
        for i in range(len(x)):

            scalars.append(pysgpp.FloatScalarVector(x[i]))

        scalars=pysgpp.FloatScalarVectorVector(scalars)
        fullgrid_operation.setParameters(scalars)

        return fullgrid_operation.eval([level for i in range(self.d)]).getValue()

    def getLevelManager(self):
        return self.LevelManager(self.operation_linear,self.level,self.d)

    class LevelManager:
        def __init__(self,operation,level,d):
            self.operation=operation
            self.level=level
            self.d=d


        #returns points with fullgrid
        def numGridPoints(self):
            return self.operation.numPoints([self.level for i in range(self.d)])
        def getAllGridPoints(self):
            return  self.operation.getGridPoints([self.level for i in range(self.d)])


class LinearFullgrid:
    def __init__(self, function, dim):
        self.func = pysgpp.multiFunc(function)
        self.grad = []
        self.level = 0
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(getgradkfunc(self.func, i)))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_linear = self.operation_linear.getFullGridEvaluator()

    def eval_full(self, level, x, fullgrid_operation):
        scalars = []
        for i in range(len(x)):
            scalars.append(pysgpp.FloatScalarVector(x[i]))

        scalars = pysgpp.FloatScalarVectorVector(scalars)
        fullgrid_operation.setParameters(scalars)

        return fullgrid_operation.eval([level for i in range(self.d)]).getValue()

    def evaluate(self, level, x):
        sum = 0
        self.level = level
        return self.eval_full(level,x,self.operation_linear)

    def getLevelManager(self):
        return self.LevelManager(self.operation_linear, self.level, self.d)

    class LevelManager:
        def __init__(self, operation, level, d):
            self.operation = operation
            self.level = level
            self.d = d

        # returns points with fullgrid
        def numGridPoints(self):
            return self.operation.numPoints([self.level for i in range(self.d)])

        def getAllGridPoints(self):
            return self.operation.getGridPoints([self.level for i in range(self.d)])



#######################################################################################################


# returns the plotgridpoints and the points wrapped as DataVector
def generate1DGrid(samples):
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, samples)
    X_eval = [pysgpp.DataVector([i]) for i in X]

    return X, X_eval


def calculate_grid_y_values(X, operation, level):
    X_eval = [pysgpp.DataVector([i]) for i in X]
    Y = [operation.evaluate(level, X_eval[i]) for i in range(len(X_eval))]

    z = operation.getLevelManager().getGridPointMatrix()

    gridpoints = [z.get(0, x) for x in range(z.getNcols())]

    return Y, gridpoints


def plot_1d_grid(X, Y, text=" ", gridpoints=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X, Y, label=text)

    if (gridpoints is not None):
        ax.plot(gridpoints, [0 for x in gridpoints], '.', c="black")

    # legend position

    return ax


def plot2DGrid(n_samples, level, operation, title=""):
    # x=0 defined 0 in combigrid
    epsilon = 10 ** -12
    X = np.linspace(0, 1, n_samples)
    Y = np.linspace(0, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = operation.evaluate(
                level, pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, z, cmap=plt.get_cmap("jet"), )
    ax.set_xlabel('$x_1$', fontsize=15)
    ax.set_ylabel('$x_2$', fontsize=15)
    ax.set_zlabel('$\~{f}$', fontsize=15)

    fig.suptitle(title)


def plot2DGrid_with_tangents(n_samples, level, operation, title=""):
    # x=0 defined 0 in combigrid
    epsilon = 10 ** -12
    X = np.linspace(0, 1, n_samples)
    Y = np.linspace(0, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = operation.evaluate(
                level, pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, z, cmap=plt.get_cmap("jet"), )

    ax.set_xlabel('$x_1$', fontsize=15)
    ax.set_ylabel('$x_2$', fontsize=15)
    ax.set_zlabel('$\~{f}$', fontsize=15)

    ### draw tangents
    operation_wrap = operationwrapper(operation, level)
    x0, x1, y = calc_gradient_tangent_gridpoint(operation.getLevelManager().getAllGridPoints(),
                                                operation_wrap, 2)



    x0=[]
    x1=[]
    y=[]
    gridpoints=operation.getLevelManager().getAllGridPoints()
    print(len(gridpoints))
    for i in range(len(gridpoints)):

        x0.append(gridpoints[i][0])
        x1.append(gridpoints[i][1])
        y.append(operation_wrap(gridpoints[i]))


    ax.scatter(x0,x1,y, c="black")

    #for i in range(len(x0)):
    #    ax.plot(x0[i], x1[i], y[i], c="gray")

    fig.suptitle(title)


def plot2D_comparison_function(n_samples, func, title=""):
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, n_samples)
    Y = np.linspace(0 + epsilon, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = func(pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.get_cmap("jet"))

    ax.set_xlabel('$x_1$', fontsize=12,labelpad=10)
    ax.set_ylabel('$x_2$', fontsize=12,labelpad=10)


    fig.suptitle(title)


def plot2DContour(n_samples, level, operation):
    epsilon = 10 ** -12
    X = np.linspace(0, 1, n_samples)
    Y = np.linspace(0, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = operation.evaluate(
                level, pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # ax.imshow(z, cmap="jet", interpolation='bilinear', extent=[1, 0, 1, 0])
    # CS = ax.contour(z, extent=[0, 1, 1, 0])
    # ax.clabel(CS, inline=1, fontsize=10)

    # grid generation
    z = operation.getLevelManager().getAllGridPoints()
    X = [z[x][0] for x in range(len(z))]
    Y = [z[y][1] for y in range(len(z))]

    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_position(('data', 1))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_position(('data', 1))

    points = ax.plot(X, Y, '.', c="black", markersize=10)


def calc_error_ilevels(operation, targetfunc, dim, maxlevel):
    nr_gridpoints = []
    errors = []
    for l in range(maxlevel):
        operation_wrap = operationwrapper(operation, l)
        gridpterror = estimatel2Error(15000, dim, operation_wrap, targetfunc)
        errors.append(gridpterror)
        nr_gridpoints.append(operation.getLevelManager().numGridPoints())
    return nr_gridpoints, errors


def plot_log_error(errors, gridpoints, text):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(errors.shape[1]):
        ax.plot(gridpoints[:, i], errors[:, i], marker='o', label=text[i])

    plt.legend()
    ax.set_xlabel("Anzahl Gitterpunkte")
    ax.set_ylabel("$L^2-Fehler$")
    ax.set_xscale("log")
    ax.set_yscale("log")


def example_1D_realfunction():
    function = f1D
    X, X_eval = generate1DGrid(500)
    Y_compare = [function(X_eval[i]) for i in range(len(X_eval))]
    plot_1d_grid(X, Y_compare, "function", )


def example_1D_psi(l):
    func = pysgpp.multiFunc(f1D)
    d = 1
    level = l
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
        d, func)
    X, _, = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    plot_1d_grid(X, Y, "$\psi$", gridpoints)



def example_1D_linear(l):
    func = pysgpp.multiFunc(f1D)
    d = 1
    level = l
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        d, func)
    X, _, = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    plot_1d_grid(X, Y, "$\linear$", gridpoints)


def example_1D_zeta(l):
    func = pysgpp.multiFunc(getgradkfunc(f1D, 0))
    d = 1
    level = l
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
        d, func)
    X, _, = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    plot_1d_grid(X, Y, "$\psi$", gridpoints)


def example_combicombigrid_1D(l):
    x = pysgpp.DataVector([0.01])
    x1 = pysgpp.DataVector([0.5])
    x2 = pysgpp.DataVector([0.99])

    d = 1
    level = l

    f1D_grad = getgradkfunc(f1D, 0)
    f1D_grad_wrap = pysgpp.multiFunc(f1D_grad)
    f1D_wrap = pysgpp.multiFunc(f1D)
    operation = CombiCombigrid1d_hermite(f1D, f1D_grad)
    op_wrap = operationwrapper(operation, level)

    grad_operation = getgradkfunc(op_wrap, 0)

    print(str(x), grad_operation(x), "ableitung")
    print(str(x1), grad_operation(x1), "ableitung")
    print(str(x1), f1D_grad(x1), "ableitung real")
    print(str(x2), grad_operation(x2), "ableitung")

    # plot the grid
    X, X_eval = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    print(
        str(operation.operation_zeta.numGridPoints()) + " grid points in zeta/psi")
    plot_1d_grid(X, Y, "combi_combigrid_1D", gridpoints)

    # error calculation

    error = estimatel2Error(
        7000, 1, operationwrapper(operation, level), f1D)
    print("the error is " + str(error))


def example_2D_linear(l, func_standard):
    d = 2
    level = l
    n_samples = 100

    operation = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        d, func_standard)

    operation_wrap = operationwrapper(operation, level)
    plot2DGrid_with_tangents(n_samples, level, operation, "linear")

    print(operation.evaluate(level, pysgpp.DataVector([0, 1])))
    print(operation.evaluate(level, pysgpp.DataVector([0.000000001, 1])))

    error = estimatel2Error(20000, 2, operation_wrap, func_standard)
    # calculate error
    print("the estimtaded l2 error for lineargrid: " + str(error))


def example_2D_psi():
    d = 2
    level = 2
    n_samples = 100
    func = pysgpp.OptRosenbrockObjective(d)
    func_standard = getfuncwrapper(func)
    func_combiwrap = pysgpp.multiFunc(func_standard)

    operation = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
        d, func_combiwrap)

    operation_wrap = operationwrapper(operation, level)
    plot2DGrid(n_samples, level, operation, func_combiwrap)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)
    # calculate error
    print("the estimtaded l2 error for psigrid: " + str(error))


def example_2D_comparison_function(func_standard,title):
    d = 2

    n_samples = 100

    plot2D_comparison_function(n_samples, func_standard, title)



def example_combicombigrid_2D_linear(l, func_standard):
    """

    Args:
        l:
        func_standard: function that is already in multifunc standard format
    """
    d = 2
    level = l
    n_samples = 50

    operation = CombiCombigridLinear(func_standard, 2)
    operation_wrap = operationwrapper(operation, level)

    plot2DGrid_with_tangents(n_samples, level, operation, "de Baar & Harding")
    plot2DContour(n_samples, level, operation)

    # derivatives


    grad_0 = getgradkfunc(operation_wrap, 0)
    grad_1 = getgradkfunc(operation_wrap, 1)

    x1 = pysgpp.DataVector([0.5, 0.5])
    # print(grad_0(x1))
    # print(grad_1(x1))


    error = estimatel2Error(20000, 2, operation_wrap, func_standard)

    # calculate error
    print("the estimtaded l2 error for combicombilinear: " + str(error))

    gridpterror = calc_error_on_gridpts(operation.getLevelManager().getAllGridPoints(),
                                        operation_wrap, func_standard)
    grad_error = calc_error_gradient(operation.getLevelManager().getAllGridPoints(), operation_wrap,
                                     func_standard, d)
    mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
                                              operation_wrap,
                                              func_standard, d)

    print("estimtaded l2 error for combicombihermite: " + str(error))
    print("1st gradient errors:" + str(grad_error))
    print("mixed gradient error:" + str(mixgrad_error))


def example_combicombigrid_2D_hermite(l, func_standard):
    d = 2
    level = l
    n_samples = 50

    operation = CombiCombigridHermite(func_standard, 2)
    #operation = HierachGridBSpline(2, 3, func_standard)
    operation_wrap = operationwrapper(operation, level)

    # operation.operation_zeta for mixed
    plot2DGrid_with_tangents(n_samples, level, operation, "combicombi_hermite")
    plot2DContour(n_samples, level, operation)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)
    grad_error = calc_error_gradient(operation.getLevelManager().getAllGridPoints(), operation_wrap,
                                     func_standard, d)
    mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
                                              operation_wrap,
                                              func_standard, d)
    # calculate error
    print("estimtaded l2 error for combicombihermite: " + str(error))
    print("1st gradient errors:" + str(grad_error))
    print("mixed gradient error:" + str(mixgrad_error))



def example_linearfullgrid(l, func_standard):
    d = 2
    level = l
    n_samples = 50

    operation = CombiFullGridHermite(func_standard,2)

    operation_wrap = operationwrapper(operation, level)

    # operation.operation_zeta for mixed
    plot2DGrid(n_samples, level, operation, "combicombi_hermite")
    plot2DContour(n_samples, level, operation)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)
    grad_error = calc_error_gradient(operation.getLevelManager().getAllGridPoints(),operation_wrap,
                                  func_standard, d)

    mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
                                             operation_wrap,
                                              func_standard, d)

    print("estimtaded l2 error for linearfullgrid: " + str(error))
    print("1st gradient errors:" + str(grad_error))
    print("mixed gradient error:" + str(mixgrad_error))


def example_plot_error(func_standard, dim,maxlevel=5):
    operation_count = 7
    operation = CombiCombigridLinear(func_standard, dim)
    operation2 = CombiCombigridHermite(func_standard, dim)
    operation3 = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        dim, func_standard)

    operation4 = CombiCombigrid2dHermite_without_mixed(func_standard)
    operation5=CombiFullGridHermite(func_standard, dim)
    operation6=CombiFullGridLinear(func_standard, dim)
    operation7=LinearFullgrid(func_standard,dim)
    # operation4= HierachGridBSpline(dim,3,func_standard)


    

    nr_gridpoints, errors = calc_error_ilevels(operation, func_standard, dim, maxlevel)
    X = np.zeros((len(errors), operation_count))
    Y = np.zeros((len(errors), operation_count), dtype=np.float64)

    X[:, 0] = nr_gridpoints
    Y[:, 0] = errors

    nr_gridpoints, errors = calc_error_ilevels(operation2, func_standard, dim, maxlevel)

    X[:, 1] = nr_gridpoints
    Y[:, 1] = errors

    nr_gridpoints, errors = calc_error_ilevels(operation3, func_standard, dim, maxlevel)

    X[:, 2] = nr_gridpoints
    Y[:, 2] = errors

    nr_gridpoints, errors = calc_error_ilevels(operation4, func_standard, dim, maxlevel)

    X[:, 3] = nr_gridpoints
    Y[:, 3] = errors

    nr_gridpoints, errors = calc_error_ilevels(operation5, func_standard, dim, maxlevel)

    X[:, 4] = nr_gridpoints
    Y[:, 4] = errors

    nr_gridpoints, errors = calc_error_ilevels(operation6, func_standard, dim, maxlevel)

    X[:, 5] = nr_gridpoints
    Y[:, 5] = errors

    nr_gridpoints, errors = calc_error_ilevels(operation7, func_standard, dim, maxlevel)

    X[:, 6] = nr_gridpoints
    Y[:, 6] = errors

    # nr_gridpoints, errors = calc_error_ilevels(operation4, func_standard, dim, maxlevel)

    # X[:, 3] = nr_gridpoints
    # Y[:, 3] = errors



    plot_log_error(Y, X, ["Gitter (de Baar und Harding)", "Gitter (Hermite)",
                          "lineares Gitter", "Gitter (Hermite) ohne gemischte Abl.",
                          "HermiteFullgrid","CombiLinearFullgrid","linearFullgrid"])


def testf(x):
    return x[0] * x[1]


# example_1D_realfunction()
#example_1D_psi(2)
#example_1D_linear(2)
# example_1D_zeta(2)
# example_combicombigrid_1D(2)


dim = 2

func = pysgpp.OptBraninObjective()
func_wrap = getfuncwrapper(func)

func_standard = pysgpp.multiFunc(func_wrap)

#example_2D_comparison_function(func_standard,"")
# example_2D_psi()
#example_2D_linear(2,func_standard)
#example_combicombigrid_2D_linear(2, func_standard)  # with "contourplot"
# print()

#example_combicombigrid_2D_hermite(3, func_standard)

example_plot_error(func_standard, dim,maxlevel=7)


#example_linearfullgrid(2,func_standard)
plt.show()

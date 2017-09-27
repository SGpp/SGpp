import pysgpp
import FunctionClasses as fctClass
import matplotlib.pyplot as plt
import numpy as np
import random
import sympy as sp

from mpl_toolkits.mplot3d import Axes3D

from itertools import chain, combinations


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


def estimatel2ErrorGradients(n, dim, gridOpEval, func_container):
    """

    Returns:
        List: the l2 error in the gradients as a list in powerset order
    """

    targetFunc = func_container.getFunction()
    mixed_directions = list(powerset([i for i in range(dim)]))

    mixed_grad = []
    opEval_mixgrad = []

    # do the derivatives for operation, function and save resulting functions in the 2 arrays
    for i in range(1, len(mixed_directions)):

        # mixed_grad_temp = targetFunc
        opEval_mixgrad_temp = gridOpEval
        for j in mixed_directions[i]:
            # mixed_grad_temp = getgradkfunc(mixed_grad_temp, j)
            opEval_mixgrad_temp = getgradkfunc(opEval_mixgrad_temp, j)
        mixed_grad_temp = func_container.getGradient(mixed_directions[i])
        mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))
        opEval_mixgrad.append(opEval_mixgrad_temp)

    errors = []
    for i in range(len(mixed_grad)):
        errors.append(estimatel2Error(n, dim, opEval_mixgrad[i], mixed_grad[i]))

    return errors


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
    mixed_directions = list(powerset([i for i in range(dim)]))

    mixed_grad = []
    opEval_mixgrad = []

    for i in range(1, len(mixed_directions)):

        mixed_grad_temp = targetFunc
        opEval_mixgrad_temp = gridOpEval
        for j in mixed_directions[i]:
            mixed_grad_temp = getgradkfunc(mixed_grad_temp, j)
            opEval_mixgrad_temp = getgradkfunc(opEval_mixgrad_temp, j)
        mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))
        opEval_mixgrad.append(opEval_mixgrad_temp)

    errors = []
    for i in range(len(mixed_grad)):
        errors.append(calc_error_on_gridpts(gridpoints, opEval_mixgrad[i], mixed_grad[i]))

    return errors


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

    if x[k] <= h:
        gy = (function(displace(x, k, h)) - function(x)) / (h)

    elif x[k] >= 1 - h:

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
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.d = dim

        self.grad = []
        for i in range(self.d):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)

        self.operation_zeta = []
        for i in range(self.d):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                    self.d, i, self.grad[i]))

    def evaluate(self, level, x):
        sum = 0
        sum += self.operation_psi.evaluate(level, x)

        for op in self.operation_zeta:
            sum += op.evaluate(level, x)

        return sum

    def getLevelManager(self):
        return self.operation_psi.getLevelManager()


class CombiCombigridHermite:
    def __init__(self, func_container, dim):
        self.d = dim
        self.func_container = func_container

        self.func = pysgpp.multiFunc(self.func_container.getFunction())
        self.grad = []

        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)

        # create all combinations for mixed gradients
        self.mixed_directions = list(powerset([i for i in range(self.d)]))
        self.mixed_grad = []
        self.operationzeta_psi = []

        for i in range(1, len(self.mixed_directions)):
            mixed_grad_temp = self.func_container.getGradient(self.mixed_directions[i])
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))
        for i in range(len(self.mixed_grad)):
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


class BaseLinearFullgrid():
    def __init__(self, dim, func_container, full=False):
        self.dim = dim

        self.func = function_eval_undisplaced(func_container.getFunction())
        self.level = -1
        self.full = full

    def linearGrid(self, level):
        grid = pysgpp.Grid.createLinearBoundaryGrid(self.dim)
        self.gridStorage = grid.getStorage()

        if (self.full == False):
            grid.getGenerator().regular(level)
        else:
            grid.getGenerator().full(level)
        alpha = pysgpp.DataVector(self.gridStorage.getSize())
        alpha.setAll(0.0)
        for i in range(self.gridStorage.getSize()):
            gp = self.gridStorage.getPoint(i)
            coordinates = pysgpp.DataVector(self.dim)
            gp.getStandardCoordinates(coordinates)
            alpha[i] = self.func.evalUndisplaced(coordinates)

        pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha)

        self.alpha = alpha
        return grid

    def evaluate(self, level, x):
        if self.level != level:
            self.level = level
            self.grid = self.linearGrid(self.level)
            self.opEvalLinear = pysgpp.createOperationEvalNaive(self.grid)

        return self.opEvalLinear.eval(self.alpha, x)

    def gridOperationEval(self, opEval, alpha):
        def eval(x):
            return opEval.eval(alpha, x)

        return eval

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


class HierachGridBSpline:
    # create a bsplineGrid with calculated hierarchized coefficients

    def __init__(self, dim, degree, func_collection, full=False):
        self.dim = dim
        self.degree = degree
        self.func = function_eval_undisplaced(func_collection.getFunction())
        self.level = -1
        self.full = full

    def createBsplineGrid(self, level):

        grid = pysgpp.Grid.createBsplineBoundaryGrid(self.dim, self.degree)
        self.gridStorage = grid.getStorage()
        print "dimensionality:         {}".format(self.gridStorage.getDimension())

        if (self.full == False):
            grid.getGenerator().regular(level)
        else:
            grid.getGenerator().full(level)

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
        if self.level != level:
            self.level = level
            self.grid, self.coeffs = self.createBsplineGrid(self.level)
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


class CombiCombigridLinear:
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

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
class FullGrid:
    def eval_full(self, level, x, fullgrid_operation):
        scalars = []
        for i in range(len(x)):
            scalars.append(pysgpp.FloatScalarVector(x[i]))

        scalars = pysgpp.FloatScalarVectorVector(scalars)
        fullgrid_operation.setParameters(scalars)

        return fullgrid_operation.eval([level for i in range(self.d)]).getValue()

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


class CombiFullGridHermite(FullGrid):
    def __init__(self, func_collection, dim):
        self.d = dim
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        self.func_collection = func_collection

        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)
        self.operation_psi = self.operation_psi.getFullGridEvaluator()

        # create all combinations for the mixed gradients
        self.mixed_directions = list(powerset([i for i in range(self.d)]))
        self.mixed_grad = []
        self.operationzeta_psi = []

        for i in range(1, len(self.mixed_directions)):
            mixed_grad_temp = self.func_collection.getGradient(self.mixed_directions[i])
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))

        for i in range(len(self.mixed_grad)):
            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    self.mixed_directions[
                                                                                        i + 1],
                                                                                    self.mixed_grad[
                                                                                        i]))
        # lazy full conversion to full grid operations

        for i in range(len(self.operationzeta_psi)):
            self.operationzeta_psi[i] = self.operationzeta_psi[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level

        sum += self.eval_full(level, x, self.operation_psi)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += self.eval_full(level, x, op)
        return sum

    def getLevelManager(self):
        return self.LevelManager(self.operation_psi, self.level, self.d)


class CombiFullGridHermite_withoutmixed(FullGrid):
    def __init__(self, func_collection, dim):
        self.d = dim
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        self.func_collection = func_collection

        # operations
        self.operation_psi = \
            pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
                self.d, self.func)
        self.operation_psi = self.operation_psi.getFullGridEvaluator()

        # create all combinations for the mixed gradients

        self.mixed_grad = []
        self.operationzeta_psi = []

        for i in range(self.d):
            mixed_grad_temp = self.func_collection.getGradient([i])
            self.mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))

        for i in range(self.d):
            self.operationzeta_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaInterpolation(self.d,
                                                                                    [i],
                                                                                    self.mixed_grad[
                                                                                        i]))
        # lazy full conversion to full grid operations

        for i in range(len(self.operationzeta_psi)):
            self.operationzeta_psi[i] = self.operationzeta_psi[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level

        sum += self.eval_full(level, x, self.operation_psi)
        # for op in self.operation_zeta_psi:
        #    sum += op.evaluate(level, x)


        for op in self.operationzeta_psi:
            sum += self.eval_full(level, x, op)
        return sum

    def getLevelManager(self):
        return self.LevelManager(self.operation_psi, self.level, self.d)


class CombiFullGridLinear(FullGrid):
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_linear = self.operation_linear.getFullGridEvaluator()

        self.operation_psi = []
        self.operation_zeta = []
        for i in range(dim):
            self.operation_psi.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryPsiLinearInterpolation(self.d,
                                                                                         i,
                                                                                         self.func))
        # lazy redefine of fullgrid
        for i in range(dim):
            self.operation_psi[i] = self.operation_psi[i].getFullGridEvaluator()

        for i in range(dim):
            self.operation_zeta.append(
                pysgpp.CombigridOperation.createExpUniformBoundaryZetaLinearInterpolation(
                    self.d, i, self.grad[i]))

        # lazy redefine of fullgrid
        for i in range(dim):
            self.operation_zeta[i] = self.operation_zeta[i].getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level

        for operation in self.operation_psi:
            sum += self.eval_full(level, x, operation)
        for operation in self.operation_zeta:
            sum += self.eval_full(level, x, operation)
        sum -= self.eval_full(level, x, self.operation_linear) * (self.d - 1)
        return sum

    def getLevelManager(self):
        return self.LevelManager(self.operation_linear, self.level, self.d)


class LinearFullgrid(FullGrid):
    def __init__(self, func_collection, dim):
        self.func = pysgpp.multiFunc(func_collection.getFunction())
        self.grad = []
        self.level = 0
        for i in range(dim):
            self.grad.append(pysgpp.multiFunc(func_collection.getGradient([i])))

        self.d = dim
        self.operation_linear = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
            self.d, self.func)

        self.operation_linear = self.operation_linear.getFullGridEvaluator()

    def evaluate(self, level, x):
        sum = 0
        self.level = level
        return self.eval_full(level, x, self.operation_linear)

    def getLevelManager(self):
        return self.LevelManager(self.operation_linear, self.level, self.d)


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


def plot2DGrid_operation(n_samples, level, operation, title=""):
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

            # if the derivative should be plotted
            # z[i,j]=getgradkfunc(operationwrapper(operation,level),0)(pysgpp.DataVector([X[i, j],
            #                                                                        Y[i, j]]))

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

    x0 = []
    x1 = []
    y = []
    gridpoints = operation.getLevelManager().getAllGridPoints()
    print(len(gridpoints))
    for i in range(len(gridpoints)):
        x0.append(gridpoints[i][0])
        x1.append(gridpoints[i][1])
        y.append(operation_wrap(gridpoints[i]))

    ax.scatter(x0, x1, y, c="black")

    # for i in range(len(x0)):
    #    ax.plot(x0[i], x1[i], y[i], c="gray")

    fig.suptitle(title)


def plot2D_comparison_function(n_samples, func, title=""):
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, n_samples)
    Y = np.linspace(0 + epsilon, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in xrange(n_samples):
        for j in range(n_samples):
            z[i, j] = func(pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.get_cmap("jet"))

    ax.set_xlabel('$x_1$', fontsize=12, labelpad=10)
    ax.set_ylabel('$x_2$', fontsize=12, labelpad=10)

    fig.suptitle(title)


####
def plot2D_comparison_function_random(n_samples, func, title=""):
    X = []
    for i in range(n_samples):
        X.append(random.random())

    Y = []
    for i in range(n_samples):
        Y.append(random.random())

    X.append(0)
    X.append(1)
    Y.append(0)
    Y.append(1)
    X.sort()
    Y.sort()

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((len(X), len(Y)))
    for i in xrange(len(X)):
        for j in range(len(Y)):
            z[i, j] = func(pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, z, cstride=1, rstride=1, cmap=plt.get_cmap("jet"))

    ax.set_xlabel('$x_1$', fontsize=12, labelpad=10)
    ax.set_ylabel('$x_2$', fontsize=12, labelpad=10)

    fig.suptitle(title)


def scatterplot2D_error(n_samples, func, title=None):
    X = []
    Y = []
    Z = []


    n_samples=(n_samples)**2

    for _ in range(n_samples):
        point = pysgpp.DataVector(2)
        for d in range(2):
            point[d] = random.random()
        # print(gridOpEval(point)-targetFunc.evalUndisplaced(point))
        X.append(point[0])
        Y.append(point[1])
        Z.append(func(point))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.title=title
    ax.scatter(X, Y, Z, c=Z,cmap="jet")


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
    for l in range(1, maxlevel):
        operation_wrap = operationwrapper(operation, l)
        gridpterror = estimatel2Error(10000, dim, operation_wrap, targetfunc)
        errors.append(gridpterror)
        nr_gridpoints.append(operation.getLevelManager().numGridPoints())
    return nr_gridpoints, errors


def calc_error_ilevels_grad(operation, func_container, dim, maxlevel, grad_index_list):
    nr_gridpoints = []
    errors = []

    for l in range(1, maxlevel):
        operation_wrap = operationwrapper(operation, l)
        for k in grad_index_list:
            operation_wrap = getgradkfunc(operation_wrap, k)
        targetfunc = func_container.getGradient(grad_index_list)

        gridpterror = estimatel2Error(10000, dim, operation_wrap, targetfunc)
        errors.append(gridpterror)
        nr_gridpoints.append(operation.getLevelManager().numGridPoints())
    return nr_gridpoints, errors


def plot_log_error(errors, gridpoints, text, x_axis_type="nr_gridpoints",title=""):
    fig = plt.figure()
    plt.title(title)
    ax = fig.add_subplot(111)
    for i in range(errors.shape[1]):
        ax.plot(gridpoints[:, i], errors[:, i], marker='o', label=text[i])

    plt.legend()
    if (x_axis_type == "nr_gridpoints"):
        ax.set_xlabel("Anzahl Gitterpunkte")
    else:
        ax.set_xlabel("Level")
    ax.set_ylabel("$L^2-Fehler$")
    if (x_axis_type == "nr_gridpoints"):
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
    plot2DGrid_operation(n_samples, level, operation, func_combiwrap)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)
    # calculate error
    print("the estimtaded l2 error for psigrid: " + str(error))


def example_2D_comparison_function(func_standard, title):
    d = 2

    n_samples = 100

    X = np.linspace(0, 1, n_samples)
    Y = np.linspace(0, 1, n_samples)

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


def example_combicombigrid_2D_hermite(l, func_collection):
    d = 2
    level = l
    n_samples = 50
    func_standard = func_collection.getFunction()
    operation = CombiCombigridLinear(func_collection, 2)
    # operation = HierachGridBSpline(2, 3, func_standard)
    operation_wrap = operationwrapper(operation, level)

    # operation.operation_zeta for mixed
    plot2DGrid_operation(n_samples, level, operation, "combicombi_hermite")
    # plot2DContour(n_samples, level, operation)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)

    # grad_error = calc_error_gradient(operation.getLevelManager().getAllGridPoints(),
    #                                 operation_wrap,
    #                                 func_standard, d)
    # mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
    #                                          operation_wrap,
    #                                          func_standard, d)
    # calculate error
    print("estimtaded l2 error for combicombihermite: " + str(error))
    # print("1st gradient errors:" + str(grad_error))
    # print("mixed gradient error:" + str(mixgrad_error))


def geterrorfunc(func1, func2):
    def errorfunc(x):
        return abs(func1(x) - func2(x))

    return errorfunc


def example_error_picewise(l, func_container):
    d = 2
    level = 5
    n_samples = 50

    operation_list = []

    func_standard = func_container.getFunction()
    func_standard = pysgpp.multiFunc(func_standard)

    operation_list.append(CombiCombigridLinear(func_container, 2))
    operation_list.append(CombiCombigridHermite(func_container, 2))
    operation_list.append(pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        2, func_standard))
    operation_list.append(CombiFullGridHermite(func_container, 2))
    operation_list.append(CombiFullGridLinear(func_container, 2))
    operation_list.append(BaseLinearFullgrid(2, func_container, True))

    text = ["Combilinear", "combiHermite", "linear", "CombiHermiteFull", "CombiLinearFull",
            "BaseLinearFull"]

    i = 0
    for operation in operation_list:
        operation_wrap = operationwrapper(operation, level)

        # operation.operation_zeta for mixed
        plot2D_comparison_function_random(n_samples, geterrorfunc(operation_wrap, func_standard),
                                          text[i])


        # plot2DContour(n_samples, level, operation)


        #error = estimatel2Error(10000, 2, operation_wrap, func_standard)
        #errorl2_gradients = estimatel2ErrorGradients(10000, 2, operation_wrap, func_container)

        # erroron_gripoint = calc_error_on_gridpts(operation.getLevelManager().getAllGridPoints(),
        #                                        operation_wrap, func_standard)
        # mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
        #                                          operation_wrap,
        #                                          func_standard, d)

        #print("estimtaded l2 error for " + text[i] + ": " + str(error))
        # print("error on gridpts: " + str(erroron_gripoint))
        # print("mixed gradient error:" + str(mixgrad_error))
        #print("l2 errors from gradients: " + str(errorl2_gradients))

        print("-----------------------------------------------------------------")
        i = i + 1


def example_plot_error(func_collection, dim, maxlevel=5, x_axis_type="nr_gridpoints",title=""):
    """

      Args:
          func_standard:
          dim:
          maxlevel:
          x_axis_type: either "nr_gridpoints" for the number of gridpoints or "level" to compare
          each operation by level
      """

    # operation4= HierachGridBSpline(dim,3,func_standard)


    func_standard = pysgpp.multiFunc(func_collection.getFunction())

    operations = []

    operations.append(CombiCombigridLinear(func_collection, dim))
    operations.append(CombiCombigridHermite(func_collection, dim))
    operations.append(pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        dim, func_standard))
    #operations.append(CombiCombigrid2dHermite_without_mixed(func_collection, dim))
    #operations.append((BaseLinearFullgrid(dim, func_collection, True)))

    operations.append(LinearFullgrid(func_collection, dim))
    operations.append(CombiFullGridHermite(func_collection, dim))
    operations.append(CombiFullGridLinear(func_collection, dim))
    # operations.append(HierachGridBSpline(dim, 3, func_collection, True))
    #operations.append(HierachGridBSpline(dim, 3, func_collection, False))

    func_standard = func_collection.getFunction()
    labels = []
    labels.append("Gitter (de Baar und Harding)")
    labels.append("Gitter (Hermite)")
    labels.append("lineares Gitter")
    #labels.append("Gitter (Hermite) ohne gemischte Abl.")

    #labels.append("BaseLinearFullgrid")
    labels.append("linearFullgrid")
    labels.append("HermiteFullgrid")
    labels.append("CombiLinearFullgrid")
    # labels.append("BSpline-Grid-Full")
    #labels.append("BSpline-Grid")

    X = np.zeros((maxlevel - 1, len(operations)))
    Y = np.zeros((maxlevel - 1, len(operations)), dtype=np.float64)

    for i in range(len(operations)):
        nr_gridpoints, errors = calc_error_ilevels(operations[i], func_standard, dim, maxlevel)

        if (x_axis_type == "level"):
            x_values = [j for j in range(1, maxlevel)]

        else:
            x_values = nr_gridpoints
        print(i, nr_gridpoints)

        X[:, i] = x_values
        Y[:, i] = errors

    plot_log_error(Y, X, labels, x_axis_type,title=title)


def example_plot_error_gradients(func_collection, dim,grad_index_list, maxlevel=5,
                                 x_axis_type="nr_gridpoints",
                                 title=""):
    """

      Args:
          func_standard:
          dim:
          maxlevel:
          x_axis_type: either "nr_gridpoints" for the number of gridpoints or "level" to compare
          each operation by level
      """

    # operation4= HierachGridBSpline(dim,3,func_standard)


    func_standard = pysgpp.multiFunc(func_collection.getFunction())

    operations = []

    operations.append(CombiCombigridLinear(func_collection, dim))
    operations.append(CombiCombigridHermite(func_collection, dim))
    operations.append(pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        dim, func_standard))
    operations.append(CombiCombigrid2dHermite_without_mixed(func_collection, dim))
    #operations.append((BaseLinearFullgrid(dim, func_collection, True)))

    operations.append(LinearFullgrid(func_collection, dim))
    operations.append(CombiFullGridHermite(func_collection, dim))
    operations.append(CombiFullGridLinear(func_collection, dim))
    operations.append(CombiFullGridHermite_withoutmixed(func_collection,dim))
    #operations.append(HierachGridBSpline(dim, 3, func_collection, True))
    #operations.append(HierachGridBSpline(dim, 3, func_collection, False))

    func_standard = func_collection.getFunction()
    labels = []
    labels.append("Gitter (de Baar und Harding)")
    labels.append("Gitter (Hermite)")
    labels.append("lineares Gitter")
    labels.append("Gitter (Hermite) ohne gemischte Abl.")

    #labels.append("BaseLinearFullgrid")
    labels.append("linearFullgrid")
    labels.append("HermiteFullgrid")
    labels.append("CombiLinearFullgrid")
    labels.append("CombiFullGridHermitewithoutmixed")
    #labels.append("BSpline-Grid-Full")
    #labels.append("BSpline-Grid")

    X = np.zeros((maxlevel-1 , len(operations)))
    Y = np.zeros((maxlevel-1 , len(operations)), dtype=np.float64)

    for i in range(len(operations)):
        nr_gridpoints, errors = calc_error_ilevels_grad(operations[i], func_container, dim,
                                                        maxlevel,grad_index_list)

        if (x_axis_type == "level"):
            x_values = [j for j in range(1, maxlevel)]

        else:
            x_values = nr_gridpoints
        print(i, nr_gridpoints)

        X[:, i] = x_values
        Y[:, i] = errors

    plot_log_error(Y, X, labels, x_axis_type,title=title)


def example_multidim(func_standard, dim, maxlevel=5):
    operation_count = 3
    operation = CombiCombigridLinear(func_standard, dim)
    operation2 = CombiCombigridHermite(func_standard, dim)
    operation3 = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        dim, func_standard)

    operation4 = CombiFullGridHermite(func_standard, dim)
    operation5 = CombiFullGridLinear(func_standard, dim)
    operation6 = LinearFullgrid(func_standard, dim)
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

    #   nr_gridpoints, errors = calc_error_ilevels(operation4, func_standard, dim, maxlevel)
    #
    #   X[:, 3] = nr_gridpoints
    #    Y[:, 3] = errors

    #   nr_gridpoints, errors = calc_error_ilevels(operation5, func_standard, dim, maxlevel)
    #
    #  X[:, 4] = nr_gridpoints
    #   Y[:, 4] = errors

    #    nr_gridpoints, errors = calc_error_ilevels(operation6, func_standard, dim, maxlevel)

    # X[:, 5] = nr_gridpoints
    # Y[:, 5] = errors

    # nr_gridpoints, errors = calc_error_ilevels(operation4, func_standard, dim, maxlevel)

    # X[:, 3] = nr_gridpoints
    # Y[:, 3] = errors



    plot_log_error(Y, X, ["Gitter (de Baar und Harding)", "Gitter (Hermite)",
                          "lineares Gitter", "HermiteFullgrid", "CombiLinearFullgrid",
                          "linearFullgrid"])


def testf(x):
    return x[0] * x[1]


def BraninSymbolic():
    x0, x1, x2 = sp.symbols("x0 x1 x2")

    x1_tmp = 15.0 * x0 - 5.0;

    x2_tmp = 15.0 * x1

    tmp = x2_tmp - 5.1 * x1_tmp * x1_tmp / (4.0 * sp.pi * sp.pi) + 5.0 * x1_tmp / sp.pi - 6.0;

    return tmp * tmp + 10.0 * (1.0 - 1.0 / (8.0 * sp.pi)) * sp.cos(x1_tmp) + 10.0;

def testfSymbolic_2d():
    x0, x1, x2 = sp.symbols("x0 x1 x2")

    x0_temp=5.0*x0
    x1_temp=5.0*x1
    return x0**2*sp.cos(x1_temp)**2*x1+(x0**2+sp.cos(x1_temp))**2

def testfSymbolic_4d():
    x0, x1, x2,x3 = sp.symbols("x0 x1 x2 x3")

    x0_temp=5.0*x0
    x1_temp=5.0*x1
    x2_temp=4.0*x2
    x3_temp=3.0*x3
    return x0**2*sp.cos(x1_temp)**2*x1*x3_temp**4+(x0**2+sp.cos(x1_temp))**2*x2_temp

def test2DSymbolic():
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    return sp.cos(x1*2)*sp.sin(0.5-x0*4)**2

def get_symbolictest():
    x0, x1, x2 = sp.symbols("x0 x1 x2")
    return (x0 * x1 + 1) ** 2 * sp.sin(x0)


# example_1D_realfunction()
# example_1D_psi(2)
# example_1D_linear(2)
# example_1D_zeta(2)
# example_combicombigrid_1D(2)


dim = 2

func = pysgpp.OptBraninObjective()
func_wrap = getfuncwrapper(func)

func_standard = pysgpp.multiFunc(func_wrap)

# example_2D_psi()
# example_2D_linear(2,func_standard)
# example_combicombigrid_2D_linear(2, func_standard)  # with "contourplot"
# print()

testclass = fctClass.funcGradientCollection(func_standard, 2)

func_container = fctClass.funcGradientCollectionSymbolic(test2DSymbolic(), 2)

#example_2D_comparison_function(func_container.getFunction(), "")

# example_combicombigrid_2D_hermite(3, func_container)
#example_plot_error(func_container, dim, maxlevel=8, x_axis_type="level")
#example_plot_error_gradients(func_container, dim,[0], maxlevel=6, x_axis_type="level",
  #                           title="x0_grad")
#example_plot_error_gradients(func_container, dim,[1], maxlevel=6, x_axis_type="level",
#                            title="x1_grad")
#example_plot_error_gradients(func_container, dim,[0,1], maxlevel=6, x_axis_type="level",
#                           title="x0_x1_grad")


#example_error_picewise(4, func_container)



dim = 4



func_container = fctClass.funcGradientCollectionSymbolic(testfSymbolic_4d(), dim)

#func = pysgpp.OptRosenbrockObjective(dim)
#func_wrap = getfuncwrapper(func)

#func_standard = pysgpp.multiFunc(func_wrap)


#example_plot_error(func_container, dim, maxlevel=5, x_axis_type="level")
#example_plot_error_gradients(func_container, dim,[0], maxlevel=6, x_axis_type="level",
#                            title="x0_grad")
example_plot_error_gradients(func_container, dim,[0,3], maxlevel=5, x_axis_type="level",
                             title="x0x3_grad")
#example_plot_error_gradients(func_container, dim,[0,1,2,3], maxlevel=5, x_axis_type="level",
#                             title="x0x1x2x3_grad")



plt.show()

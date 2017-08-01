import pysgpp

import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D


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


def generateGridpoints_1D(level):
    gridpoints = []

    if (level is 0):
        return [0.5]

    i = 0
    while i * 2 ** (-level) <= 1:
        gridpoints.append(i * 2 ** (-level))
        i += 1

    return gridpoints


def generateGridpoints_2D(level):

    gp1D=generateGridpoints_1D(level)
    X,Y=np.meshgrid(gp1D,gp1D)
    return X,Y



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
    h = 1e-09

    gy = (function(displace(x, k, h)) - function(displace(x, k, -h))) / (2 * h)
    return gy


# parabola between [0,1]
def f1D(x):
    return -4 * ((x[0] - 0.5) ** 2) + 1

def f1D_grad(x):
    return 4-8*x[0]


def f2D(x):
    return 4 * x[0] * (1 - x[0]) * 4 * x[1] * (1 - x[1])



def f2D_test(x):
    return x[0]**2+x[1]**2

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


# not ready
class CombiCombigrid2dHermite:
    def __init__(self, function):
        self.func = pysgpp.multiFunc(function)
        self.grad_0 = pysgpp.multiFunc(getgradkfunc(self.func, 0))
        self.grad_1 = pysgpp.multiFunc(getgradkfunc(self.func, 1))

        self.d = 2
        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)
        self.operation_psi_zeta_0 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, 0, self.grad_0)
        self.operation_psi_zeta_1 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, 0, self.grad_1)
        self.operation_zeta = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
            self.d, self.grad_func)

    def evaluate(self, level, x):
        sum = 0
        sum += self.operation_psi.evaluate(level, x)
        sum += self.operation_psi_zeta_0.evaluate(level,
                                                  x) + self.operation_psi_zeta_1(
            level, x)
        sum += self.operation_zeta(level, x)
        return sum


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


# returns the plotgridpoints and the points wrapped as DataVector
def generate1DGrid(samples):
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, samples)
    X_eval = [pysgpp.DataVector([i]) for i in X]

    return X, X_eval


def calculate_grid_y_values(X, operation, level):
    X_eval = [pysgpp.DataVector([i]) for i in X]
    Y = [operation.evaluate(level, X_eval[i]) for i in range(len(X_eval))]

    gridpoints = generateGridpoints_1D(level)

    return Y, gridpoints


def plot_1d_grid(X, Y, text=" ", gridpoints=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X, Y, label=text)

    if (gridpoints is not None):
        ax.plot(gridpoints, [0 for x in gridpoints], '.', c="black")

    # legend position

    return ax



def plot2DGrid(n_samples, level, operation,title=""):
    # x=0 defined 0 in combigrid
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, n_samples)
    Y = np.linspace(0 + epsilon, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = operation.evaluate(
                level, pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.get_cmap("jet"),)

    fig.suptitle(title)


def plot2D_comparison_function(n_samples,func,title=""):
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

    fig.suptitle(title)


def plot2DContour(n_samples, level, operation):
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, n_samples)
    Y = np.linspace(0 + epsilon, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = operation.evaluate(
                level, pysgpp.DataVector([X[i, j], Y[i, j]]))


    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(z,cmap="jet",interpolation='bilinear',extent=[1,0,1,0])
    CS=ax.contour(z,extent=[0,1,1,0])
    ax.clabel(CS, inline=1, fontsize=10)


    X,Y=generateGridpoints_2D(level)



    points=ax.plot(X,Y,'.',c="black",markersize=10)


def example_1D_realfunction():
    function=f1D
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






def example_2D_linear(l,func_standard):
    d = 2
    level = l
    n_samples = 100


    operation = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        d, func_standard)

    operation_wrap=operationwrapper(operation,level)
    plot2DGrid(n_samples, level, operation, "linear")

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

    operation_wrap=operationwrapper(operation,level)
    plot2DGrid(n_samples, level, operation, func_combiwrap)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)
    # calculate error
    print("the estimtaded l2 error for psigrid: " + str(error))


def example_2D_comparison_function(func_standard):
    d = 2

    n_samples = 100

    plot2D_comparison_function(n_samples, func_standard,"real function")


def example_combicombigrid_2D_linear(l, func_standard):
    """

    Args:
        l:
        func_standard: function that is already in multifunc standard format
    """
    d = 2
    level = l
    n_samples =50





    operation = CombiCombigrid2dLinear(func_standard)
    operation_wrap=operationwrapper(operation,level)


    plot2DGrid(n_samples,level,operation,"combicombi_linear")
    plot2DContour(n_samples, level, operation)

    # derivatives


    grad_0=getgradkfunc(operation_wrap,0)
    grad_1=getgradkfunc(operation_wrap, 1)

    x1 = pysgpp.DataVector([0.5,0.5])
    print(grad_0(x1))
    print(grad_1(x1))


    error=estimatel2Error(20000,2,operation_wrap,func_standard)
    # calculate error
    print("the estimtaded l2 error for combicombilinear: "+str(error))


def example_combicombigrid_2D_hermite(l):
    d = 2
    level = l
    n_samples =50
    func = pysgpp.OptRosenbrockObjective(d)

    func_standard = getfuncwrapper(func)


    operation = CombiCombigrid2dHermite(func_standard)
    operation_wrap=operationwrapper(operation,level)


    x = pysgpp.DataVector([0.4, 0.5])

    plot2DGrid(n_samples,level,operation,"combicombi_linear")
    plot2DContour(n_samples, level, operation)


    error=estimatel2Error(10000,2,operation_wrap,func_standard)
    # calculate error
    print("the estimtaded l2 error for combicombilinear: "+str(error))



# example_1D_realfunction()
# example_1D_psi(3)
# example_combicombigrid_1D(3)


func = pysgpp.OptBubbleWrapObjective(2)
func_standard=pysgpp.multiFunc(func)

#getfuncwrapper(func)

example_2D_comparison_function(func_standard)
#example_2D_psi()
example_2D_linear(3,func_standard)
example_combicombigrid_2D_linear(3,func_standard)
#example_combicombigrid_2D_hermite(3)


plt.show()

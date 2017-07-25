import pysgpp
import math
import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D


def operationwrapper(combi_operation, level):
    def operation(x):
        return combi_operation.evaluate(level, x)
    return operation

# for functions from sgpp::optimization
def getfuncwrapper(func, dim):
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
        sum += (gridOpEval(point) - targetFunc(point))**2
    sum = sum / n
    sum = np.sqrt(sum)
    return sum







def displace(x, i, h):

    x_temp = pysgpp.DataVector(x)
    x_temp[i] = x[i] + h
    return x_temp

#get the derivatives for point x
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
    return -4 * ((x[0] - 0.5)**2) + 1
    
def f2D(x):
    return 4 * x[0] * (1 - x[0]) * 4 * x[1] * (1 - x[1])



#combi_combigrids

class combi_combigrid_1D:
    def __init__(self, function, grad_function):
        self.func = pysgpp.multiFunc(function)
        self.grad_func = pysgpp.multiFunc(grad_function)
        self.d = 1
        self.operation_psi = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)
        self.operation_zeta = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
            self.d, self.grad_func)

    def evaluate(self, level, x):
        return self.operation_psi.evaluate(level, x) + 1 * self.operation_zeta.evaluate(level, x)


# not ready
class combi_combigrid_2D_hermite:
    def __init__(self, function, grad_function):
        self.func = pysgpp.multiFunc(function)
        self.grad_func = pysgpp.multiFunc(grad_function)
        self.d = 2
        self.operation_psi = ysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, self.func)
        self.operation_zeta_0 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, 0, self.func)
        self.operation_zeta_1 = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
            self.d, 0, self.func)
        self.operation_zeta = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
            self.d, self.grad_func)

    def evaluate(self, level, x):
        return self.operation_psi.evaluate(level, x) + self.operation_zeta.evaluate(level, x)



class combi_combigrid_2D_linear:
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





def generate1DGrid(samples, level, operation):
    epsilon = 10**-12
    X = np.linspace(0 + epsilon, 1, samples)
    X_eval = [pysgpp.DataVector([i]) for i in X]
    #print(operation.evaluate(level, X_eval[0]))

    Y = [operation.evaluate(level, X_eval[i]) for i in range(len(X_eval))]

    return X, X_eval, Y

def plot1DGrid(X, Y):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X,  Y, label='$\psi$-evaluation')


def plot2DGrid(n_samples, level, operation, func=None):


    # x=0 defined 0 in combigrid
    epsilon = 10**-12
    X = np.linspace(0 + epsilon, 1, n_samples)
    Y = np.linspace(0 + epsilon, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = operation.evaluate(
                level, pysgpp.DataVector([X[i, j], Y[i, j]]))
            #z[i, j] = func(pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, z, cmap=plt.get_cmap("jet"))
    


def example_2D_psi():

    d = 2
    level = 4
    n_samples = 100
    testfunc = pysgpp.OptAckleyObjective(d)
    f2D = getfuncwrapper(testfunc, d)
    func = pysgpp.multiFunc(f2D)

    operation = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
        d, func)
    plot2DGrid(n_samples, level, operation, func)


def example_1D_psi():
    func = pysgpp.multiFunc(f1D)
    d = 1
    level = 3
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
        d, func)
    X,_, Y = generate1DGrid(500, level, operation)
    plot1DGrid(X, Y)


def example_combicombigrid_1D():
    x = pysgpp.DataVector([0.5])
    x1 = pysgpp.DataVector([0.0])
    x2 = pysgpp.DataVector([1.0])
    k = 0
    d = 1
    level = 2

    f1D_grad = getgradkfunc(f1D, 0)
    f1D_grad_wrap = pysgpp.multiFunc(f1D_grad)
    f1D_wrap = pysgpp.multiFunc(f1D)
    operation = combi_combigrid_1D(f1D, f1D_grad)
    # operation=pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(1,f1D_wrap)

    print(f1D_grad(x), "ableitung (0,5)")
    print(f1D_grad(x1), "ableitung (0)")
    print(f1D_grad(x2), "ableitung (1)")

    # print(operation.operation_zeta.evaluate(0,x))

    # plot the grid
    X, X_eval, Y = generate1DGrid(500, level, operation)
    print(operation.operation_zeta.numGridPoints(), "grid points in zeta/psi")
    plot1DGrid(X, Y)

    # plot real values
    Y_compare = [f1D(X_eval[i])for i in range(len(X_eval))]
    plot1DGrid(X, Y_compare)

    # error calculation

    print("the error is", estimatel2Error(
        7000, 1, operationwrapper(operation, level), f1D))


def example_combicombigrid_2D_linear():

    d = 2
    level = 4
    n_samples = 100
    testfunc = pysgpp.OptRosenbrockObjective(d)
    f2D = getfuncwrapper(testfunc, d)

    operation = combi_combigrid_2D_linear(f2D)

    plot2DGrid(n_samples, level, operation, f2D)
    x = pysgpp.DataVector([0.4, 0.5])

    print(operation.evaluate(0, x))


example_1D_psi()
#example_2D_psi()
example_combicombigrid_1D()
#example_combicombigrid_2D_linear()
plt.show()

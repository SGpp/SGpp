import pysgpp
import FunctionClasses as fctClass
import matplotlib.pyplot as plt
import numpy as np
import GridCollection as gC
import random
import sympy as sp
from mpl_toolkits.mplot3d import Axes3D


def operationwrapper(combi_operation, level):
    def operation(x):
        return combi_operation.evaluate(level, x)

    return operation



def plot_1d_grid(X, Y, text=" ", gridpoints=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X, Y, label=text)

    if (gridpoints is not None):
        ax.plot(gridpoints, [0 for x in gridpoints], '.', c="black")

    # legend position

    return ax


def plot2DGrid_operation(n_samples, level, operation, title="",show=True):
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
    if(show):
        plt.show()


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
    x0, x1, y = fctClass.calc_gradient_tangent_gridpoint(operation.getLevelManager(

    ).getAllGridPoints(),
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


def plot2D_comparison_function(n_samples, func, title="",show=True):

    X = np.linspace(0 , 1, n_samples)
    Y = np.linspace(0 , 1, n_samples)

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

    if (show == True):
        plt.show()


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

    ax.imshow(z, cmap="jet", interpolation='bilinear', extent=[1, 0, 1, 0])
    CS = ax.contour(z, extent=[0, 1, 1, 0])
    ax.clabel(CS, inline=1, fontsize=10)

    # grid generation
    z = operation.getLevelManager().getAllGridPoints()
    X = [z[x][0] for x in range(len(z))]
    Y = [z[y][1] for y in range(len(z))]

    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_position(('data', 1))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_position(('data', 1))

    points = ax.plot(X, Y, '.', c="black", markersize=10)


# gets a function and not a operation
def plot2DContour_func(n_samples, func,operation,title=""):

    X = np.linspace(0, 1, n_samples)
    Y = np.linspace(0, 1, n_samples)

    X, Y = np.meshgrid(X, Y)

    z = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            z[i, j] = func(pysgpp.DataVector([X[i, j], Y[i, j]]))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(z, cmap="jet", interpolation='bilinear', extent=[0, 1, 0, 1])
    CS = ax.contour(z, extent=[0, 1, 1, 0])
    ax.clabel(CS, inline=1, fontsize=10)

    # grid generation
    z = operation.getLevelManager().getAllGridPoints()
    X = [z[x][0] for x in range(len(z))]
    Y = [z[y][1] for y in range(len(z))]

    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_position(('data', 1))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_position(('data', 1))

    points = ax.plot(X, Y, '.', c="white", markersize=10,)
    ax.set_title(title)
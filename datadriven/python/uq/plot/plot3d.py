from pysgpp.extensions.datadriven.uq.operations import evalSGFunction
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from pysgpp import DataVector
from scipy.stats import gaussian_kde

import matplotlib.pyplot as plt
import numpy as np


def plotDensity3d(U, n=50):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    xlim = [0, 1]  # U.getBounds()[0]
    ylim = [0, 1]  # U.getBounds()[1]
    X = np.linspace(xlim[0], xlim[1], n)
    Y = np.linspace(ylim[0], ylim[1], n)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros(n * n).reshape(n, n)

    for i in xrange(len(X)):
        for j, (x, y) in enumerate(zip(X[i], Y[i])):
            Z[i, j] = U.pdf([x, y])

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    # ax.set_zlim(0, 2)

    fig.colorbar(surf, shrink=0.5, aspect=5)

    return fig, ax


def plotSG3d(grid, alpha, n=50, f=lambda x: x):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(0, 1, n)
    Y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros(n * n).reshape(n, n)
    for i in xrange(len(X)):
        for j, (x, y) in enumerate(zip(X[i], Y[i])):
            Z[i, j] = f(evalSGFunction(grid, alpha, DataVector([x, y])))

    # get grid points
    gs = grid.getStorage()
    gps = np.zeros([gs.size(), 2])
    p = DataVector(2)
    for i in xrange(gs.size()):
        gs.get(i).getCoords(p)
        gps[i, :] = p.array()

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.scatter(gps[:, 0], gps[:, 1], np.zeros(gs.size()))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # ax.set_zlim(0, 2)

    fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig, ax, Z


def plotFunction3d(f, xlim=[0, 1], ylim=[0, 1], n=50):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(xlim[0], xlim[1], n)
    Y = np.linspace(ylim[0], ylim[1], n)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros(n * n).reshape(n, n)

    for i in xrange(len(X)):
        for j, (x, y) in enumerate(zip(X[i], Y[i])):
            Z[i, j] = f([x, y])

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    # ax.set_zlim(0, 2)

    fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig, ax, Z


def plotSGNodal3d(grid, alpha):
    gs = grid.getStorage()

    A = np.ndarray([gs.size(), 3])

    p = DataVector(2)
    for i in xrange(gs.size()):
        gs.get(i).getCoords(p)
        A[i, 0] = p[0]
        A[i, 1] = p[1]
        A[i, 2] = evalSGFunction(grid, alpha, p)

    return plotNodal3d(A), A


def plotNodal3d(A):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(A[:, 0], A[:, 1], A[:, 2])

    return fig, ax


def plotKDE3d(values):
    kde = gaussian_kde(values)
    density = kde(values)

    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
    x, y = values
    ax.scatter(x, y, density, color='red')
    return fig, ax

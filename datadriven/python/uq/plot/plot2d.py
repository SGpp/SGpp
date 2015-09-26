from pysgpp.extensions.datadriven.uq.operations import evalSGFunction
from pysgpp import DataVector, DataMatrix

import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunctionMulti


def plotDensity2d(U, n=50, addContour=True):
    xlim, ylim = [0, 1], [0, 1]  # U.getBounds()

    x = np.linspace(xlim[0], xlim[1], n)
    y = np.linspace(ylim[0], ylim[1], n)
    X, Y = np.meshgrid(x, y)
    Z = np.ones((n, n))

    for i in xrange(len(X)):
        for j, (xi, yi) in enumerate(zip(X[i], Y[i])):
            Z[i, j] = U.pdf([xi, 1 - yi])

    # np.savetxt('density2d.csv', z.reshape(n * n, 3), delimiter=' ')

    plt.imshow(Z, interpolation='bicubic', aspect='auto',
               extent=[xlim[0], xlim[1], ylim[0], ylim[1]])

    plt.jet()
    plt.colorbar()

    if addContour:
        cs = plt.contour(X, 1 - Y, Z, colors='black')
        plt.clabel(cs, inline=1, fontsize=18)

def plotSGDE2d(U, n=100):
    gs = U.grid.getStorage()
    x = [0.0] * gs.size()
    y = [0.0] * gs.size()

    for i in xrange(gs.size()):
        x[i] = gs.get(i).getCoord(0)
        y[i] = gs.get(i).getCoord(1)

    neg_x = []
    neg_y = []
    neg_z = []
    xlim, ylim = U.getBounds()
    for xi in np.linspace(xlim[0], ylim[1], n):
        for yi in np.linspace(ylim[0], ylim[1], n):
            value = U.pdf([yi, 1 - xi])
            if value < 0 and abs(value) > 1e-14:
                neg_x.append(yi)
                neg_y.append(1 - xi)
                neg_z.append(U.pdf([yi, 1 - xi]))

    # plot image of density
    plotDensity2d(U)

    # plot grid points
#     plt.plot(x, y, linestyle=' ', marker='o', color="white")

    # plot data points
#     plt.plot(U.trainData.array()[:, 0], U.trainData.array()[:, 1],
#              linestyle=' ', marker='x', color='orange')

    # plot negative areas
    if len(neg_z) > 0:
        plt.plot(neg_x, neg_y, linestyle=' ', marker='o', color='red')
        plt.title("N=%i, [%g, %g]" % (U.grid.getStorage().size(), min(neg_z), max(neg_z)))
    else:
        plt.title("N=%i" % U.grid.getStorage().size())
    plt.xlim(0, 1)
    plt.ylim(0, 1)


def plotFunction2d(f, addContour=True, n=101):
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y)
    Z = np.ones(n * n).reshape(n, n)

    print "-" * 60
    for i in xrange(len(X)):
        for j, (xi, yi) in enumerate(zip(X[i], Y[i])):
            Z[i, j] = f(xi, yi)
            if yi == 0.6:
                print xi, yi, f(xi, yi)

    plt.imshow(Z, interpolation='bilinear', extent=(0, 1, 0, 1))

    plt.jet()
    plt.colorbar()

    if addContour:
        cs = plt.contour(X, 1 - Y, Z, colors='black')
        plt.clabel(cs, inline=1, fontsize=18)

    return


def plotSG2d(grid, alpha, addContour=True, n=100):
    gs = grid.getStorage()

    gpxp = []
    gpyp = []

    gpxn = []
    gpyn = []

    for i in xrange(gs.size()):
        if alpha[i] > 0:
            gpxp.append(gs.get(i).getCoord(0))
            gpyp.append(gs.get(i).getCoord(1))
        else:
            gpxn.append(gs.get(i).getCoord(0))
            gpyn.append(gs.get(i).getCoord(1))

    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y)
    Z = np.ones(n * n).reshape(n, n)

    neg_x = []
    neg_y = []
    neg_z = []

    A = DataMatrix(n * n, 2)
    p = DataVector(2)
    # do vectorized evaluation
    k = 0
    for i in xrange(len(X)):
        for j, (xi, yi) in enumerate(zip(X[i], Y[i])):
            p[0] = xi
            p[1] = 1 - yi
            A.setRow(k, p)
            k += 1

    res = evalSGFunctionMulti(grid, alpha, A)

    k = 0
    for i in xrange(len(X)):
        for j, (xi, yi) in enumerate(zip(X[i], Y[i])):
            Z[i, j] = res[k]
            if Z[i, j] < 0 and abs(Z[i, j]) > 1e-13:
                neg_x.append(xi)
                neg_y.append(1 - yi)
                neg_z.append(res[k])
            k += 1

    plt.imshow(Z, interpolation='bilinear', extent=(0, 1, 0, 1))

    if len(neg_z) > 0:
        plt.plot(neg_x, neg_y, linestyle=' ', marker='o', color='red')
        plt.title("[%g, %g]" % (min(neg_z), max(neg_z)))

    # plot surpluses
    plt.plot(gpxp, gpyp, "^ ", color="white")
    plt.plot(gpxn, gpyn, "v ", color="red")

    plt.jet()
    plt.colorbar()

    if addContour:
        cs = plt.contour(X, 1 - Y, Z, colors='black')
        plt.clabel(cs, inline=1, fontsize=18)

    return res


def plotGrid2d(grid, alpha=None):
    gs = grid.getStorage()
    gps = {'p': np.zeros([0, 2]),
           'n': np.zeros([0, 2])}
    p = DataVector(2)
    for i in xrange(gs.size()):
        gs.get(i).getCoords(p)
        if alpha is None or alpha[i] >= 0:
            gps['p'] = np.vstack((gps['p'], p.array()))
        else:
            gps['n'] = np.vstack((gps['n'], p.array()))

    # plot the grid points
    plt.plot(gps['p'][:, 0], gps['p'][:, 1], "^ ", color='red')
    plt.plot(gps['n'][:, 0], gps['n'][:, 1], "v ", color='red')
    plt.xlim(0, 1)
    plt.ylim(0, 1)


def plotSamples2d(samples):
    X = np.zeros(len(samples))
    Y = np.zeros(len(samples))
    for i, sample in enumerate(samples):
        X[i], Y[i] = sample.getActiveProbabilistic()
    plt.plot(X, Y, linestyle=' ', marker='o')

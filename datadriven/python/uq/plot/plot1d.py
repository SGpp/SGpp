import numpy as np
import matplotlib.pyplot as plt
from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.operations import evalSGFunction
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import dehierarchize
from pysgpp.extensions.datadriven.uq.plot.colors import loadColorSequence


def plotFunction1d(f, n=1000, xlim=[0, 1], **kws):
    x = np.linspace(xlim[0], xlim[1], n, endpoint=True)
    y = [f(np.array([xi])) for xi in x]
    plt.plot(x, y, **kws)

def plotCDF1d(U, n=1000, *args, **kws):
    bounds = U.getBounds()
    if len(bounds) == 1:
        bounds = bounds[0]
    x = np.linspace(bounds[0], bounds[1], n)
    y = [U.cdf(xi) for xi in x]

    plt.plot(x, y, *args, **kws)

def plotDensity1d(U, n=1000, *args, **kws):
    bounds = U.getBounds()
    if len(bounds) == 1:
        bounds = bounds[0]
    x = np.linspace(bounds[0], bounds[1], n)
    y = [U.pdf(xi) for xi in x]

    plt.plot(x, y, *args, **kws)

def plotSGDE1d(U, n=1000):
    x = np.linspace(0, 1, n, endpoint=True)
    y = [U.pdf(xi) for xi in x]
    plt.plot(x, y)


def plotGrid1d(grid, f=lambda x: x):
    gs = grid.getStorage()
    x = np.zeros(gs.getSize())
    nodalValues = np.zeros(gs.getSize())
    for i in xrange(gs.getSize()):
        x[i] = f(gs.getCoordinate(gs.getPoint(i), 0))

    plt.scatter(x, nodalValues, marker="o")


def plotNodal1d(A):
    plt.plot(A[:, 0], A[:, 1], " ", marker="o")

def plotSGNodal1d(grid, alpha):
    gs = grid.getStorage()
    A = np.ndarray([gs.getSize(), 2])

    p = DataVector(2)
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        A[i, 0] = p[0]
        A[i, 1] = evalSGFunction(grid, alpha, p.array())

    return plotNodal1d(A), A

def plotSG1d(grid, alpha, n=1000, f=lambda x: x, show_grid_points=False,
             **kws):
    x = np.linspace(0, 1, n)
    y = [f(evalSGFunction(grid, alpha, np.array([xi])))
         for xi in x]

    plt.plot(x, y, **kws)

    if show_grid_points:
        gs = grid.getStorage()
        for i in xrange(gs.getSize()):
            x = gs.getCoordinate(gs.getPoint(i), 0)
            plt.scatter(x, 0, color="black")


def plotCDF(p, buckets):
    fig = plt.figure()
    # p, buckets = U.cdf()
    # plt.bar(buckets[:-1] + np.diff(buckets) / 2, p[1:])
    plt.plot(buckets, p, ' ', marker='o')
    plt.ylim(0, 1)
    return fig
    # x = np.linspace(buckets.min(), buckets.max(), 1000)
    # y = [self.sol_pdf.cdf(xi) for xi in x]
    # plt.plot(x, y)


def plotPDF(p, buckets):
    fig = plt.figure()
    plt.plot(buckets, p, ' ', marker='o')
    plt.ylim(0, 1)
    return fig


def plotSurplusLevelWise(data, maxLevel):
    fig = plt.figure()
    for level, surpluses in data.items():
        plt.plot([level] * len(surpluses), surpluses, ' ', marker='o')
    plt.xlim(np.min(data.keys()) - 1, maxLevel + 1)
    return fig


def plotSobolIndices(sobolIndices, ts=None, legend=False,
                     names=None):
    fig = plt.figure()
    plots = []
    colors = loadColorSequence(len(sobolIndices))

    if legend and names is None:
        raise Exception("plotSobolIndices - attribute names is not set")

    if ts is None:
        y0 = 0
        for i in xrange(len(sobolIndices)):
            myplot = plt.bar([0], [sobolIndices[i]], 1, bottom=[y0], color=colors[i])
            y0 += sobolIndices[i]
            plots = [myplot] + plots

        if legend:
            plt.xticks([0.5], ('sobol indices',))
            plt.ylim(0, 1)
            plt.xlim(-0.2, 2)
            plt.legend(plots,
                       [r"$S_{%s}$ = %.3f" % (name, value)
                        for (name, value) in zip(names[::-1],
                                                 sobolIndices[::-1])])

    else:
        y0 = np.zeros(sobolIndices.shape[0])
        for i in xrange(sobolIndices.shape[1]):
            y1 = y0 + sobolIndices[:, i]
            myplot, = plt.plot(ts, y1, lw=4)
            plt.fill_between(ts, y0, y1, color=colors[i % len(colors)], alpha=.5)
            y0 = y1

            plots = [myplot] + plots

        if legend:
            plt.xlim(min(ts), max(ts) + (max(ts) - min(ts)) / 6.)
            plt.legend(plots,
                       [r"$S_{%s}$" % (name,) for name in names[::-1]],
                       loc='upper right')
    return fig

from pysgpp import DataVector
from pysgpp.extensions.datadriven.uq.operations import evalSGFunction

import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import dehierarchize


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


def plotNodal1d(grid, alpha):
    gs = grid.getStorage()
    x = np.zeros(gs.size())
    nodalValues = np.zeros(gs.size())
    for i in xrange(gs.size()):
        x[i] = gs.get(i).getCoord(0)
        nodalValues[i] = evalSGFunction(grid, alpha, DataVector([x[i]]))

    plt.plot(x, nodalValues, " ", marker="o")


def plotSG1d(grid, alpha, n=1000, f=lambda x: x, **kws):
    x = np.linspace(0, 1, n)
    y = [f(evalSGFunction(grid, alpha, np.array([xi])))
         for xi in x]

    plt.plot(x, y, **kws)


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
    return fig


def plotSobolIndices(sobolIndices, ts=None, legend=False,
                     names=None,
                     colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']):
    fig = plt.figure()
    plots = []

    if legend and names is None:
        raise Exception("plotSobolIndices - attribute names is not set")

    if ts is None:
        y0 = 0
        for i in xrange(len(sobolIndices)):
            myplot = plt.bar([0], [sobolIndices[i]], 1, color=colors[i], bottom=[y0])
            y0 += sobolIndices[i]
            plots = [myplot] + plots

        if legend:
            plt.xticks([0.5], ('sobol indices',))
            plt.ylim(0, 1)
            plt.xlim(-0.2, 2)
            plt.legend(plots,
                       [r"$S_{%s}$ = %g" % (name, value)
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

import numpy as np
import matplotlib.pyplot as plt
from colors import load_color, load_font_properties
from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.operations import evalSGFunction
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import dehierarchize
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend


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


def plotHistogram1d(samples,
                    mean_label=None,
                    mean_color=load_color(1),
                    std_label=None,
                    std_color=load_color(2),
                    median_label=None,
                    median_color=load_color(3),
                    *args, **kws):
    n, bins, patches = plt.hist(samples,
                                normed=True, color=load_color(0),
                                edgecolor="white", *args, **kws)

    if mean_label is not None:
        mean = np.mean(samples)
        ix = np.where(bins <= mean)[0]
        if len(ix) > 0:
            ix = ix[-1]
        else:
            ix = 0
        plt.vlines(mean, 0.0, n[ix], color=mean_color, label=mean_label,
                   linewidth=2)
    if std_label is not None:
        mean = np.mean(samples)
        std_err = np.std(samples)
        ix = np.where(bins <= mean - std_err)[0]
        if len(ix) > 0:
            ix = ix[-1]
        else:
            ix = 0
        plt.vlines(mean - std_err, 0.0, n[ix], color=std_color, label=std_label,
                   linewidth=2)
        ix = np.where(bins <= mean + std_err)[0][-1]
        plt.vlines(mean + std_err, 0.0, n[ix], color=std_color, linewidth=2)
    if median_label is not None:
        median = np.median(samples)
        ix = np.where(bins <= median)[0]
        if len(ix) > 0:
            ix = ix[-1]
        else:
            ix = 0
        plt.vlines(median, 0.0, n[ix], color=median_color, label=median_label,
                   linewidth=2)

def plotDensity1d(U, n=1000,
                  alpha_value=None,
                  mean_label=None,
                  facecolor=load_color(1),
                  interval_label=None,
                  *args, **kws):
    bounds = U.getBounds()
    if len(bounds) == 1:
        bounds = bounds[0]
    x = np.linspace(bounds[0], bounds[1], n)
    y = np.array([U.pdf(xi) for xi in x])

    plt.plot(x, y, *args, **kws)

    if mean_label is not None:
        plt.vlines(U.mean(), 0.0, U.pdf(U.mean()), color=facecolor, label=mean_label)

    if alpha_value is not None:
        # define label for interval plot
        if interval_label is None:
            interval_label = r"$[F(\alpha / 2), F(1 - \alpha/2)]$"
        # show interval that contains 1 - alpha
        x_min, x_max = U.ppf(alpha_value / 2.), U.ppf(1. - alpha_value / 2.)
        ixs = np.intersect1d(np.where(x_min <= x),
                             np.where(x <= x_max))
#         plt.vlines(x_min, 0.0, y[ixs.min()], color=facecolor, label=interval_label)
#         plt.vlines(x_max, 0.0, y[ixs.max()], color=facecolor)
        plt.fill_between(x[ixs],
                         y[ixs],
                         np.zeros(y[ixs].shape[0]),
                         facecolor=facecolor, alpha=0.2,
                         label=interval_label)


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
                     adjust_yaxis=True, names=None, mc_reference=None):
    fig = plt.figure()
    plots = []

    if legend and names is None:
        raise Exception("plotSobolIndices - attribute names is not set")

    lgd = None
    if ts is None:
        y0 = 0
        for i in xrange(len(sobolIndices)):
            myplot = plt.bar([0], [sobolIndices[i]], 1, bottom=[y0], color=load_color(i))
            y0 += sobolIndices[i]
            plots = [myplot] + plots

        if legend:
            plt.xticks([0.5], ('sobol indices',))
            if adjust_yaxis:
                plt.ylim(0, 1)

            plt.xlim(-0.2, 2)
            lgd = plt.legend(plots,
                             [r"$S_{%s}$ = %.3f" % (name, value)
                              for (name, value) in zip(names[::-1],
                                                       sobolIndices[::-1])],
                             prop=load_font_properties())

    else:
        y0 = np.zeros(sobolIndices.shape[0])
        offset = 1 if mc_reference is not None else 0
        for i in xrange(sobolIndices.shape[1]):
            y1 = y0 + sobolIndices[:, i]
            color = load_color(i + offset)
            myplot, = plt.plot(ts, y1, color=color, lw=2)
            plt.fill_between(ts, y0, y1, color=color, alpha=.5)
            y0 = y1

            plots = [myplot] + plots

        labels = [r"$S_{%s}$" % (",".join(name),) for name in names[::-1]]
        if mc_reference is not None:
            myplot, = plt.plot(mc_reference["ts"],
                               mc_reference["values"],
                               marker=mc_reference["marker"],
                               color=mc_reference["color"])
            plt.fill_between(mc_reference["ts"],
                             mc_reference["values"],
                             mc_reference["err"][:, 0],
                             facecolor=mc_reference["color"], alpha=0.2)
            plt.fill_between(mc_reference["ts"],
                             mc_reference["values"],
                             mc_reference["err"][:, 1],
                             facecolor=mc_reference["color"], alpha=0.2)
            labels = [mc_reference["label"]] + labels
            plots = [myplot] + plots

        if legend:
            plt.xlim(min(ts), max(ts))

            if adjust_yaxis:
                plt.ylim(0, 1)

            fig.tight_layout()
            ax = plt.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
            lgd = plt.legend(plots,
                             labels,
                             loc='upper left',
                             bbox_to_anchor=(1.02, 1),
                             borderaxespad=0,
                             prop=load_font_properties())

    return fig, lgd

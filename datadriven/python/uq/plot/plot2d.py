import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
from matplotlib.mlab import griddata

from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.operations import evalSGFunction
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunctionMulti


def addContours(xv, yv, Z,
                levels=None,
                clabels=None,
                manual_locations=None):
    if levels is not None:
        cs = plt.contour(xv, yv, Z,
                         levels=levels,
                         colors='white')
        if clabels is None:
            clabels = levels

        fmt = {}
        for level in clabels:
            if level % 1 == 0:
                fmt[level] = "%1.0f" % level
            elif level * 10 % 1 == 0:
                fmt[level] = "%1.1f" % level
            elif level * 100 % 1 == 0:
                fmt[level] = "%1.2f" % level
            else:
                fmt[level] = "%1.3f" % level

        plt.clabel(cs, clabels,
                   inline=1,
                   color="white",
                   fontsize=18,
                   fmt=fmt,
                   manual=manual_locations)

        #             # add white rectangles under clabels
        #             rect = Rectangle((10 * (24 * 60 * 60), 0.6),
#                              11 * (24 * 60 * 60), 0.8,
#                              facecolor="white",
#                              alpha=0.95, zorder=5)
#             ax = plt.gca()
#             ax.add_patch(rect)
    else:
        cs = plt.contour(xv, yv, Z, colors='white')
        plt.clabel(cs, inline=1, fmt="%1.2f", fontsize=18)


def plotTimedependentDensity2dWithRawData(xv, yv, Z, ts, us,
                                          addContour=True,
                                          color_bar_label=r'$\hat{F}(\xi_1, \xi_2)$',
                                          levels=None,
                                          clabels=None,
                                          manual_locations=None):
    # np.savetxt('density2d.csv', z.reshape(n * n, 3), delimiter=' ')
    im = plt.imshow(Z,
                    interpolation='bicubic',
                    origin="lower",
                    aspect='auto',
                    extent=[ts.min(), ts.max(), us.min(), us.max()])

    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel(color_bar_label)
    plt.clim(0, 1)

    if addContour:
        addContours(xv, yv, Z, levels, clabels, manual_locations)


def plotTimedependentDensity2d(Us, us, ts,
                               addContour=True,
                               color_bar_label=r'$\hat{F}(\xi_1, \xi_2)$',
                               levels=None,
                               clabels=None,
                               manual_locations=None):
    Z = np.ones((us.shape[0], ts.shape[0]))

    xv, yv = np.meshgrid(ts, us, sparse=False, indexing='xy')
    for i in xrange(len(ts)):
        for j in xrange(len(us)):
            Z[j, i] = Us[i](np.array([yv[j, i]]))

    # np.savetxt('density2d.csv', z.reshape(n * n, 3), delimiter=' ')
    im = plt.imshow(Z,
                    interpolation='bicubic',
                    origin="lower",
                    aspect='auto',
                    extent=[ts.min(), ts.max(), us.min(), us.max()])

    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel(color_bar_label)
    plt.clim(0, 1)

    if addContour:
        addContours(xv, yv, Z, levels, clabels, manual_locations)


def plotDensity2d(U, n=50, addContour=True,
                  color_bar_label=r'$\hat{f}(\xi_1, \xi_2)$',
                  levels=None,
                  clabels=None,
                  manual_locations=None):
    xlim, ylim = U.getBounds()

    x = np.linspace(xlim[0], xlim[1], n)
    y = np.linspace(ylim[0], ylim[1], n)
    Z = np.ones((n, n))

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            Z[j, i] = U.pdf(np.array([xv[j, i], yv[j, i]]))

    # np.savetxt('density2d.csv', z.reshape(n * n, 3), delimiter=' ')
    im = plt.imshow(Z,
                    interpolation='bicubic',
                    origin="lower",
                    aspect='auto',
                    extent=[xlim[0], xlim[1], ylim[0], ylim[1]])

    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel(color_bar_label)

    if addContour:
        addContours(xv, yv, Z, levels, clabels, manual_locations)


def plotSGDE2d(U, n=100):
    gs = U.grid.getStorage()
    x = [0.0] * gs.getSize()
    y = [0.0] * gs.getSize()

    for i in xrange(gs.getSize()):
        x[i] = gs.getCoordinate(gs.getPoint(i), 0)
        y[i] = gs.getCoordinate(gs.getPoint(i), 1)

    neg_x = []
    neg_y = []
    neg_z = []
    xlim, ylim = U.getBounds()

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            value = U.pdf([xv[j, i], yv[j, i]])

            if value < 0 and abs(value) > 1e-14:
                neg_x.append(xv[j, i])
                neg_y.append(yv[j, i])
                neg_z.append(value)

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
        plt.title("N=%i, [%g, %g]" % (U.grid.getSize(), min(neg_z), max(neg_z)))
    else:
        plt.title("N=%i" % U.grid.getSize())

    plt.xlim(0, 1)
    plt.ylim(0, 1)


def plotFunction2d(f, addContour=True, n=101,
                   xlim=[0, 1], ylim=[0, 1],
                   color_bar_label=r'$u(\xi_1, \xi_2)$',
                   levels=None,
                   clabels=None,
                   manual_locations=None):
    x = np.linspace(xlim[0], xlim[1], n)
    y = np.linspace(ylim[0], ylim[1], n)
    Z = np.ones((n, n))

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            Z[j, i] = f(np.array([xv[j, i], yv[j, i]]))

    # np.savetxt('density2d.csv', z.reshape(n * n, 3), delimiter=' ')
    im = plt.imshow(Z[::-1, :], interpolation='bicubic', aspect='auto',
                    extent=[xlim[0], xlim[1], ylim[0], ylim[1]])
    cbar = plt.colorbar(im)
    cbar.ax.set_ylabel(color_bar_label)

    if addContour:
        addContours(xv, yv, Z, levels, clabels, manual_locations)


def plotSG2d(grid, alpha, addContour=True, n=100,
             show_negative=False, show_grid_points=False,
             show_numbers=False,
             colorbarLabel=r"$\hat{f}_{\mathcal{I}}(\xi)$",
             levels=None,
             clabels=None,
             manual_locations=None):
    gs = grid.getStorage()

    gpxp = []
    gpyp = []

    gpxn = []
    gpyn = []

    gpxz = []
    gpyz = []

    numbers = []

    p = DataVector(2)
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)

        if alpha[i] > 1e-14:
            gpxp.append(p[0])
            gpyp.append(p[1])
        elif alpha[i] < -1e-14:
            gpxn.append(p[0])
            gpyn.append(p[1])
        else:
            gpxz.append(p[0])
            gpyz.append(p[1])

        numbers.append((i, p[0], p[1]))

    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    Z = np.ones(n * n).reshape(n, n)

    neg_x = []
    neg_y = []
    neg_z = []

    A = np.ndarray((n * n, 2))
    # do vectorized evaluation
    k = 0
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            A[k, :] = [xv[j, i], yv[j, i]]
            k += 1

    res = evalSGFunctionMulti(grid, alpha, A)

    k = 0
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            Z[j, i] = res[k]
            if Z[j, i] < 0 and abs(Z[j, i]) > 1e-13:
                neg_x.append(xv[j, i])
                neg_y.append(yv[j, i])
                neg_z.append(res[k])
            k += 1

    im = plt.imshow(Z[::-1, :], interpolation='bilinear', extent=(0, 1, 0, 1))

    if len(neg_z) > 0 and show_negative:
        plt.plot(neg_x, neg_y, linestyle=' ', marker='o', color='red')
        plt.title("[%g, %g]" % (min(neg_z), max(neg_z)))

    # plot surpluses
    if show_grid_points:
        plt.plot(gpxp, gpyp, "^ ", color="white")
        plt.plot(gpxn, gpyn, "v ", color="red")
        plt.plot(gpxz, gpyz, "o ", color="white")

    if show_numbers:
        for i, x, y in numbers:
            plt.text(x, y, "%i" % i, color='yellow', fontsize=12)

    cbar = plt.colorbar(im)
    cbar.set_label(colorbarLabel)

    if addContour:
        addContours(xv, yv, Z, levels, clabels, manual_locations)

    return res


def plotGrid2d(grid, alpha=None, show_numbers=True,
               xlim=(0, 1), ylim=(0, 1),
               *args, **kws):
    gs = grid.getStorage()
    gps = {'a': np.ndarray((0, 2)),
           'p': np.ndarray((0, 2)),
           'n': np.ndarray((0, 2))}

    p = DataVector(2)
    numbers = []
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        if alpha is None:
            gps['a'] = np.vstack((gps['a'], p.array()))
        else:
            if alpha[i] >= 0:
                gps['p'] = np.vstack((gps['p'], p.array()))
            else:
                gps['n'] = np.vstack((gps['n'], p.array()))

        numbers.append((i, p[0], p[1]))

    # plot the grid points
    if alpha is None:
        plt.plot(gps['a'][:, 0], gps['a'][:, 1], "o ", color='black', *args, **kws)
    else:
        plt.plot(gps['p'][:, 0], gps['p'][:, 1], "^ ", color='blue', *args, **kws)
        plt.plot(gps['n'][:, 0], gps['n'][:, 1], "v ", color='red', *args, **kws)
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])

    if show_numbers:
        for i, x, y in numbers:
            plt.text(x, y, "%i" % i, color='black', fontsize=12)


def plotSamples2d(samples):
    X = np.zeros(len(samples))
    Y = np.zeros(len(samples))
    for i, sample in enumerate(samples):
        X[i], Y[i] = sample.getActiveProbabilistic()
    plt.plot(X, Y, linestyle=' ', marker='o')


def plotScatter2d(samples, values, bounds=None,
                  color_bar_label=r"$u(x)$"):
    # define grid.
    if bounds is None:
        x_min, x_max = samples[:, 0].min(), samples[:, 0].max()
        y_min, y_max = samples[:, 1].min(), samples[:, 1].max()
    else:
        x_min, x_max = bounds[0]
        y_min, y_max = bounds[1]

    xs = np.linspace(x_min, x_max, 100)
    ys = np.linspace(y_min, y_max, 100)
    # grid the data.
    zs = griddata(samples[:, 0],
                  samples[:, 1],
                  values,
                  xs, ys,
                  interp="nn")
    # contour the gridded data, plotting dots at the randomly spaced data points.
    cs = plt.contour(xs, ys, zs, 15, linewidths=0.5, colors='k')
    cs = plt.contourf(xs, ys, zs, 15)

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel(color_bar_label)

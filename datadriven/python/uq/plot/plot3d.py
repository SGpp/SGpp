from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

from pysgpp import DataVector, createOperationEval, createOperationEvalNaive
from pysgpp.extensions.datadriven.uq.operations import evalSGFunction

from colors import load_default_color_map


def plotDensity3d(U, n=36):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    xlim = U.getBounds()[0]
    ylim = U.getBounds()[1]
    x = np.linspace(xlim[0], xlim[1], n + 1, endpoint=True)
    y = np.linspace(ylim[0], ylim[1], n + 1, endpoint=True)
    Z = np.zeros((n + 1, n + 1))

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            Z[j, i] = U.pdf(np.array([xv[j, i], yv[j, i]]))

    ax.plot_wireframe(xv, yv, Z, color="black")
    cset = ax.contour(xv, yv, Z, zdir='z', offset=np.min(Z), cmap=cm.coolwarm)
    cset = ax.contour(xv, yv, Z, zdir='x', offset=xlim[0], cmap=cm.coolwarm)
    cset = ax.contour(xv, yv, Z, zdir='y', offset=ylim[1], cmap=cm.coolwarm)

#     surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
#                            cmap=cm.coolwarm, linewidth=0, antialiased=False)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    # ax.set_zlim(0, 2)

#     fig.colorbar(surf, shrink=0.5, aspect=5)

    return fig, ax, Z


def plotGrid3d(grid, grid_points_at=0, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    # get grid points
    gs = grid.getStorage()
    gps = np.zeros([gs.getSize(), 2])
    p = DataVector(2)
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        gps[i, :] = p.array()

    ax.plot(gps[:, 0], gps[:, 1], np.ones(gps.shape[0]) * grid_points_at,
           " ", c="red", marker="o", ms=15)


def insert_labels_3d(ax,
                     xlabel=r"$\xi_1$",
                     ylabel=r"$\xi_2$",
                     zlabel=r"$u(\xi_1, \xi_2)$"):
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.xaxis.labelpad = 13
    ax.yaxis.labelpad = 13
    ax.zaxis.labelpad = 12


def plotSG3d(grid, alpha, n=36,
             f=lambda x: x, grid_points_at=0, z_min=np.Inf,
             isConsistent=True, show_grid=True,
             surface_plot=False,
             xoffset=0.0, yoffset=1.0):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = np.linspace(0, 1, n + 1, endpoint=True)
    y = np.linspace(0, 1, n + 1, endpoint=True)
    Z = np.zeros((n + 1, n + 1))

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            Z[j, i] = f(evalSGFunction(grid, alpha,
                                       np.array([xv[j, i], yv[j, i]]),
                                       isConsistent=isConsistent))

    if surface_plot:
        cmap = load_default_color_map()
        norm = colors.Normalize(vmin = Z.min(),
                                vmax = Z.max(),
                                clip = False)
        ax.plot_surface(xv, yv, Z,
                        rstride=1, cstride=1,
                        norm=norm, cmap=cmap)
    else:
        ax.plot_wireframe(xv, yv, Z, color="black")

    z_min, z_max = min(np.min(Z), z_min), np.max(Z)
    ax.set_zlim(z_min, z_max)
    if np.any(np.abs(alpha) > 1e-13):
        cset = ax.contour(xv, yv, Z, zdir='z', offset=z_min, cmap=cm.coolwarm)
        cset = ax.contour(xv, yv, Z, zdir='x', offset=xoffset, cmap=cm.coolwarm)
        cset = ax.contour(xv, yv, Z, zdir='y', offset=yoffset, cmap=cm.coolwarm)

    if show_grid:
        # get grid points
        gs = grid.getStorage()
        gps = np.zeros([gs.getSize(), 2])
        p = DataVector(2)
        for i in xrange(gs.getSize()):
            gs.getCoordinates(gs.getPoint(i), p)
            p0, p1 = p.array()
            color = "red" if alpha[i] < 0 else "blue"
            ax.plot(np.array([p0]),
                    np.array([p1]),
                    np.array([grid_points_at]),
                    " ", c=color, marker="o", ms=15)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # ax.set_zlim(0, 2)

#     fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig, ax, Z



def plotFunction3d(f, xlim=[0, 1], ylim=[0, 1], n=36,
                   z_min=np.Inf, xoffset=0.0, yoffset=1.0):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = np.linspace(xlim[0], xlim[1], n + 1, endpoint=True)
    y = np.linspace(ylim[0], ylim[1], n + 1, endpoint=True)
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    Z = np.zeros((n + 1, n + 1))

    for i in xrange(len(x)):
        for j in xrange(len(y)):
            Z[j, i] = f(np.array([xv[j, i], yv[j, i]]))

    ax.plot_wireframe(xv, yv, Z, color="black")

    z_min, z_max = min(np.min(Z), z_min), np.max(Z)
    ax.set_zlim(z_min, z_max)
    cset = ax.contour(xv, yv, Z, zdir='z', offset=z_min, cmap=cm.coolwarm)
    cset = ax.contour(xv, yv, Z, zdir='x', offset=xoffset, cmap=cm.coolwarm)
    cset = ax.contour(xv, yv, Z, zdir='y', offset=yoffset, cmap=cm.coolwarm)

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    # ax.set_zlim(0, 2)

#     fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig, ax, Z


def plotSGNodal3d(grid, alpha):
    gs = grid.getStorage()

    A = np.ndarray([gs.getSize(), 3])

    p = DataVector(2)
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        A[i, 0] = p[0]
        A[i, 1] = p[1]
        A[i, 2] = evalSGFunction(grid, alpha, p.array())

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

def plotError3d(f1, f2, xlim=[0, 1], ylim=[0, 1], n=32):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = np.linspace(xlim[0], xlim[1], n + 1, endpoint=True)
    y = np.linspace(ylim[0], ylim[1], n + 1, endpoint=True)
    Z = np.zeros((n + 1, n + 1))

    xv, yv = np.meshgrid(x, y, sparse=False, indexing='xy')
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            xi = np.array([xv[j, i], yv[j, i]])
            Z[j, i] = np.abs(f1(xi) - f2(xi))

    ax.plot_wireframe(xv, yv, Z, color="black")
    cset = ax.contour(xv, yv, Z, zdir='z', offset=np.min(Z), cmap=cm.coolwarm)
    cset = ax.contour(xv, yv, Z, zdir='x', offset=xlim[0], cmap=cm.coolwarm)
    cset = ax.contour(xv, yv, Z, zdir='y', offset=ylim[1], cmap=cm.coolwarm)

#     surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#                            linewidth=0, antialiased=False)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    # ax.set_zlim(0, 2)

#     fig.colorbar(surf, shrink=0.5, aspect=5)
    return fig, ax, Z

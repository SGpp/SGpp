# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import matplotlib.pyplot as plt

from pysgpp import Grid, DataVector
from pysgpp.extensions.datadriven.uq.operations import extend_grid
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hasBorder


def plotGrid2d(grid):
    gs = grid.getStorage()
    fig = plt.figure()
    plt.plot()
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    p = DataVector(2)
    for i in range(gs.getSize()):
        gp = gs.getPoint(i)
        gs.getCoordinates(gp, p)
        plt.plot(p[0], p[1], marker='o', color='blue')
        plt.text(p[0], p[1], "(%i, %i), (%i, %i)" % (gp.getLevel(0),
                                                     gp.getLevel(1),
                                                     gp.getIndex(0),
                                                     gp.getIndex(1)),
                 horizontalalignment='center')
    fig.show()


def plotGrid3dSlices(grid):
    gs = grid.getStorage()
    p = DataVector(3)

    d = {}

    for i in range(gs.getSize()):
        gp = gs.getPoint(i)
        gs.getCoordinates(gp, p)
        if p[2] in d:
            d[p[2]].append([p[0], p[1],
                            gp.getLevel(0), gp.getLevel(1),
                            gp.getIndex(0), gp.getIndex(1)])
        else:
            d[p[2]] = [[p[0], p[1],
                        gp.getLevel(0), gp.getLevel(1),
                        gp.getIndex(0), gp.getIndex(1)]]

    print(sum([len(dd) for dd in list(d.values())]))

    for z, items in list(d.items()):
        fig = plt.figure()
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title('z = %g %s, len=%i' % (z, "(border)" if hasBorder(grid) else "",
                                         len(items)))
        for x, y, l0, l1, i0, i1 in items:
            plt.plot(x, y, marker='o', color='blue')
            plt.text(x, y, "(%i, %i), (%i, %i)" % (l0, l1, i0, i1),
                     horizontalalignment='center')
        fig.show()


def plot_1d_2d(level=2):
    grid1d = Grid.createLinearGrid(1)
    grid1d.getGenerator().regular(level)
    grid2d2 = extend_grid(grid1d, 1)
    plotGrid2d(grid2d2)

    grid2d = Grid.createLinearGrid(2)
    grid2d.getGenerator().full(level)
    plotGrid2d(grid2d)

    plt.show()


def plot_2d_3d(level=2):
    grid2d = Grid.createLinearGrid(2)
    grid2d.getGenerator().regular(level)
    plotGrid2d(grid2d)
    grid3d = extend_grid(grid2d, 1)
    plotGrid3dSlices(grid3d)

    plt.show()


plot_1d_2d()
plot_2d_3d()

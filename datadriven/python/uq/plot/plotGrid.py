from pysgpp import createOperationEval, DataVector

import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunction
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d


def plotGrid(grid, alpha, admissibleSet, params, refined=None):
    gs = grid.getStorage()
    T = params.getJointTransformation()
    p = DataVector(2)

    x = np.ndarray((gs.getSize(), 2))
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        x[i, :] = T.unitToProbabilistic(p.array())


    a = np.ndarray((gs.getSize(), 2))
    for i, gp in enumerate(admissibleSet):
        gs.getCoordinates(gp, p)
        a[i, :] = T.unitToProbabilistic(p.array())

    r = np.ndarray((len(refined), 2))
    if refined:
        for i, gpi in enumerate(refined):
            gs.getCoordinates(gpi, p)
            r[i, :] = T.unitToProbabilistic(p.array())

    n = 50
    U = params.getIndependentJointDistribution()
    plotDensity2d(U)

    # plot grid points
    plt.plot(x[:, 0], x[:, 1], linestyle=' ', marker='o', color='g', markersize=20)  # grid
    plt.plot(a[:, 0], a[:, 1], linestyle=' ', marker='^', color='y', markersize=20)  # admissible set
    plt.plot(r[:, 0], r[:, 1], linestyle=' ', marker='v', color='r', markersize=20)  # refined points
    plt.title("size = %i" % gs.getSize())

    global myid
    # plt.savefig("out_%i.jpg" % (myid))
    # plt.close()
    myid += 1

myid = 0

from pysgpp import createOperationEval, DataVector

import numpy as np
import matplotlib.pyplot as plt


def plotGrid(grid, alpha, admissibleSet, params, refined=None):
    gs = grid.getStorage()
    x = [0.0] * gs.getSize()
    y = [0.0] * gs.getSize()

    for i in xrange(gs.getSize()):
        x[i] = gs.get(i).getCoord(0)
        y[i] = gs.get(i).getCoord(1)

    xa = [0.0] * len(admissibleSet)
    ya = [0.0] * len(admissibleSet)
    for i, gp in enumerate(admissibleSet):
        xa[i] = gp.getCoord(0)
        ya[i] = gp.getCoord(1)

    xr = []
    yr = []
    if refined:
        xr = [0.0] * len(refined)
        yr = [0.0] * len(refined)

        for i, ix in enumerate(refined):
            xr[i] = gs.get(ix).getCoord(0)
            yr[i] = gs.get(ix).getCoord(1)

    n = 50
    A = np.ones(n * n).reshape(n, n)
    B = np.ones(n * n).reshape(n, n)
    U = params.getIndependentJointDistribution()
    bounds = U.getBounds()
    opEval = createOperationEval(grid)
    for i, xi in enumerate(np.linspace(bounds[0][0], bounds[0][1], n)):
        for j, yj in enumerate(np.linspace(bounds[1][0], bounds[1][1], n)):
            A[i, j] = U.pdf([yj, 1 - xi])
            B[i, j] = opEval.eval(alpha, DataVector([xi, yj]))

    fig = plt.figure()
    plt.imshow(A, interpolation='bicubic', extent=[0,1,0,1])

    plt.jet()
    plt.colorbar()

    plt.plot(x, y, linestyle=' ', marker='o', color='g', markersize=20)     # grid
    # plt.plot(xa, ya, linestyle=' ', marker='^', color = 'y', markersize=20) # admissible set
    plt.plot(xr, yr, linestyle=' ', marker='v', color = 'r', markersize=20) # refined points
    plt.title("size = %i" % gs.getSize())
    # plt.xlim(0, 1)
    # plt.ylim(0, 1)

    # fig = plt.figure()
    # plt.imshow(B, interpolation='bicubic', extent=[0,1,0,1])

    # plt.jet()
    # plt.colorbar()
    # fig.show()
    global myid
    # plt.savefig("out_%i.jpg" % (myid))
    # plt.close()
    myid += 1

    return fig

myid = 0

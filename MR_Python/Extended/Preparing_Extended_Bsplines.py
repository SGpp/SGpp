import ipdb
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def func1(x):
    return  1


def funcx(x):
    return  x[0] 


def funcx2(x):
    return  x[0] * x[0] 


def funcx3(x):
    return  x[0] * x[0] * x[0]


def funcx4(x):
    return  x[0] * x[0] * x[0] * x[0]


def funcx5(x):
    return  x[0] * x[0] * x[0] * x[0] * x[0]


def interpolateAndError(degree, dim, level, numErrPoints, func):
    grid = pysgpp.Grid.createNakBsplineExtendedGrid(dim, degree)
    # grid = pysgpp.Grid.createNakBsplineBoundaryGrid(dim, degree)
    grid.getGenerator().regular(level)
    gridStorage = grid.getStorage()
    f_values = pysgpp.DataVector(gridStorage.getSize()) 
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        p = np.zeros(dim)
        for d in range(dim):
            p[d] = gp.getStandardCoordinate(d)
        f_values[i] = func(p)
    alpha = pysgpp.DataVector(len(f_values))
    hierSLE = pysgpp.HierarchisationSLE(grid)
    sleSolver = pysgpp.AutoSLESolver()
    if not sleSolver.solve(hierSLE, f_values, alpha):
        print "Solving failed, exiting."
        sys.exit(1)
    I = pysgpp.InterpolantScalarFunction(grid, alpha)
     
    x = pysgpp.DataVector(dim, 0) 
    err = np.zeros(numErrPoints + 1)
    points = np.zeros(numErrPoints + 1)
    for i in range(numErrPoints + 1):
        x[0] = float(i) / numErrPoints 
        points[i] = x[0]
        err[i] = func(x) - I.eval(x)
    return np.linalg.norm(err)


pysgpp.Printer.getInstance().setVerbosity(-1)
dim = 1
numErrPoints = 1000
degree = 5
for degree in [1, 3, 5]:
    print("degree {}".format(degree))
    for level in range(1, 5):
        err = np.zeros(6)
        err [0] = interpolateAndError(degree, dim, level, numErrPoints, func1)
        err [1] = interpolateAndError(degree, dim, level, numErrPoints, funcx)
        err [2] = interpolateAndError(degree, dim, level, numErrPoints, funcx2)
        err [3] = interpolateAndError(degree, dim, level, numErrPoints, funcx3)
        err [4] = interpolateAndError(degree, dim, level, numErrPoints, funcx4)
        err [5] = interpolateAndError(degree, dim, level, numErrPoints, funcx5)
        print("{}: {}".format(level, err))

# plt.figure()
# plt.plot(points, err)

#---------- plotting ----------
plt.figure()
degree = 5
level = 5
# I = [ 0, 1, 2, 3, 4, 5, 6, 7, 8]
I = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
basis = pysgpp.SNakBsplineExtendedBase(degree)
# basis = pysgpp.SNakBsplineBase(degree)
# X = np.linspace(0 - 3.0 / (2 ** level), 1 + 3.0 / (2 ** level), 1000)
X = np.linspace(0, 1, 100)
B = np.zeros((len(I), len(X)))
for j, x in enumerate(X):
    for i in range(len(I)):
        B[i, j] = basis.eval(level, I[i], x)
for i in range(len(I)):         
    plt.plot(X, B[i, :])
plt.show()

import pysgpp
import sys


def f(x, y):
    return x+y


dim = 2
# B-spline degree
p = 3
# sparse grid level
level = 3


grid = pysgpp.Grid_createOterSplineBoundaryGrid(2, 3)
gridStorage = grid.getStorage()
grid.getGenerator().regular(level)

print(f'Number of grid points: {grid.getSize()}')

# interpolation
functionValues = pysgpp.DataVector(grid.getSize())
for i in range(gridStorage.getSize()):
    gp = gridStorage.getPoint(i)
    functionValues[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))

coeffs = pysgpp.DataVector(len(functionValues))
hierSLE = pysgpp.HierarchisationSLE(grid)
sleSolver = pysgpp.Eigen()
# solve linear system
if not sleSolver.solve(hierSLE, functionValues, coeffs):
    print("Solving failed, exiting.")
    sys.exit(1)

# create interpolant
ft = pysgpp.InterpolantScalarFunction(grid, coeffs)

# evalaute interpolant
print(ft.eval(pysgpp.DataVector([0.33, 0.77])))

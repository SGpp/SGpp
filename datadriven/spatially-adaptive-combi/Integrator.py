import abc,logging
import numpy as np
import scipy as sp
import scipy.integrate
from scipy.interpolate import interpn
import math

# This is the abstract interface of an integrator that integrates a given area specified by start for function f
# using numPoints many points per dimension
class IntegratorBase(object):
    @abc.abstractmethod
    def __call__(self, f, numPoints, start, end):
        pass


# This integrator computes the trapezoidal rule for the given interval without constructing the grid explicitly
class IntegratorTrapezoidalFast(IntegratorBase):
    def __call__(self, f, numPoints, start, end):
        dim = len(start)
        length = np.empty(dim)
        offsets = np.ones(dim, dtype=np.int64)
        spacing = np.empty(dim)
        gridsize = np.int64(1)
        for i in range(dim):
            length[i] = end[i] - start[i]
            spacing[i] = float(length[i]) / float(numPoints[i] - 1)
            gridsize *= np.int64(numPoints[i])
            if i != 0:
                offsets[i] = offsets[i - 1] * int(numPoints[i - 1])
        h_prod = np.prod(spacing)
        result = 0.0
        for i in range(gridsize):
            position = np.zeros(dim)
            rest = i
            factor = 0
            for d in reversed(list(range(dim))):
                position[d] = start[d] + int(rest / offsets[d]) * spacing[d]
                if int(rest / offsets[d]) == 0 or int(rest / offsets[d]) == numPoints[d] - 1:
                    factor += 1
                rest = rest % offsets[d]
            result += f(position) * 0.5 ** factor * h_prod
        del length
        del offsets
        del spacing
        return result


# This integrator computes the trapezoidal rule for the given interval by constructing the grid explicitly
# and applying iteratively 1D trapezoidal rules
class IntegratorGenerateGridTrapezoidal(IntegratorBase):
    def __call_(self, f, numPoints, start, end):
        dim = len(start)
        length = np.zeros(dim)
        offsets = np.ones(dim, dtype=np.int64)
        spacing = np.zeros(dim)
        gridsize = np.int64(1)
        for i in range(dim):
            length[i] = end[i] - start[i]
            spacing[i] = float(length[i]) / float(numPoints[i] - 1)
            gridsize *= np.int64(numPoints[i])
            if i != 0:
                offsets[i] = offsets[i - 1] * int(numPoints[i - 1])
        startTime = time.time()
        gridValues = np.zeros(gridsize)
        for i in range(gridsize):
            position = np.zeros(dim)
            rest = i
            for d in reversed(list(range(dim))):
                position[d] = start[d] + int(rest / offsets[d]) * spacing[d]
                rest = rest % offsets[d]
            gridValues[i] = f(position)
        endTime = time.time()
        startTime = time.time()
        currentSliceSize = gridsize
        for d in reversed(list(range(dim))):
            currentSliceSize = int(currentSliceSize / int(numPoints[d]))
            for i in range(currentSliceSize):
                lineValues = np.zeros(int(numPoints[d]))
                for j in range(int(numPoints[d])):
                    lineValues[j] = gridValues[i + j * offsets[d]].copy()
                gridValues[i] = np.trapz(lineValues, dx=spacing[d])
                del lineValues
        endTime = time.time()
        result = gridValues[0].copy()
        del gridValues
        gc.collect()
        return result


# This integrator computes the integral of an arbitrary grid from the Grid class
# using the predefined interfaces and weights. The grid is not explicitly constructed.
class IntegratorArbitraryGrid(IntegratorBase):
    def __init__(self, grid):
        self.grid = grid

    def __call__(self, f, numPoints, start, end):
        dim = len(start)
        offsets = np.ones(dim, dtype=np.int64)
        gridsize = np.int64(1)
        for i in range(dim):
            gridsize *= np.int64(numPoints[i])
            if i != 0:
                offsets[i] = offsets[i - 1] * int(numPoints[i - 1])
        result = 0.0
        for i in range(gridsize):
            indexvector = np.empty(dim, dtype=int)
            rest = i
            for d in range(dim - 1, -1, -1):
                indexvector[d] = int(rest / offsets[d])
                rest = rest % offsets[d]
            result += self.integrate_point(f, indexvector)
        del offsets
        return result

    def integrate_point(self, f, indexvector):
        position = self.grid.getCoordinate(indexvector)
        return f(position) * self.grid.getWeight(indexvector)


# This integrator computes the integral of an arbitrary grid from the Grid class
# using the predefined interfaces and weights. The grid is explicitly constructed and efficiently evaluated using numpy.
class IntegratorArbitraryGridScalarProduct(IntegratorBase):
    def __init__(self, grid):
        self.grid = grid

    def __call__(self, f, numPoints, start, end):
        points, weights = self.grid.get_points_and_weights()
        f_values = [f(point) for point in points]
        return np.inner(f_values, weights)

'''
#This integrator computes the integral of an arbitrary grid from the Grid class
#using the predefined interfaces and weights. The grid is not explicitly constructed.
#GridPointsAreExcluded
class IntegratorArbitraryGridNoBoundary(IntegratorBase):
    def __init__(self,grid):
        self.grid = grid

    def __call__(self,f,numPoints,start,end):
        dim = len(start)
        offsets = np.ones(dim,dtype=np.int64)
        gridsize = np.int64(1)
        numPointsEval = list(numPoints)
        indexOffset = np.zeros(dim,dtype=np.int64)
        for i in range(dim):
            if start[i]==0 :
                numPointsEval[i] -= 1
                indexOffset[i] = 1
            if end[i] == 1:
                numPointsEval[i] -= 1
            gridsize *= np.int64(numPointsEval[i])
            if i != 0:
                offsets[i] = offsets[i-1] * int(numPointsEval[i-1])
        result = 0.0
        for i in range(gridsize):
            indexvector = np.zeros(dim,dtype=int)
            rest = i
            for d in reversed(list(range(dim))):
                indexvector[d] = int(rest / offsets[d]  + indexOffset[d])
                rest = rest % offsets[d]
            position = self.grid.getCoordinate(indexvector)
            result += f(position) * self.grid.getWeight(indexvector)
        del offsets
        return result
'''


# This is a helper method used in the single dimension method.
# It interpolates the grid to a finer grid and computes partial integrals for the subareas
def integrateVariableStartNumpyArbitraryDimAndInterpolate(f, numPoints, start, end, numberOfGridsContained):
    dim = len(start)
    length = np.zeros(dim)
    offsets = np.ones(dim, dtype=np.int64)
    extendedOffsets = np.ones(dim, dtype=np.int64)  # extended to interpolated Array

    spacing = np.zeros(dim)
    extendedSpacing = np.zeros(dim)
    gridsize = np.int64(1)
    for i in range(dim):
        length[i] = end[i] - start[i]
        spacing[i] = float(length[i]) / float(numPoints[i] - 1)
        extendedSpacing[i] = spacing[i] / numberOfGridsContained[i]
        gridsize *= np.int64(numPoints[i])
        if i != 0:
            offsets[i] = offsets[i - 1] * int(numPoints[i - 1])
            extendedOffsets[i] = offsets[i - 1] * (((int(numPoints[i - 1]) - 1) * numberOfGridsContained[i - 1]) + 1)
    startTime = time.time()

    gridValues = np.zeros(gridsize)
    for i in range(gridsize):
        position = np.zeros(dim)
        rest = i
        for d in reversed(list(range(dim))):
            position[d] = start[d] + int(rest / offsets[d]) * spacing[d]
            rest = rest % offsets[d]
        gridValues[i] = f(position)
    # number of Points of extended Array
    extendedNumPoints = [(numPoints[d] - 1) * numberOfGridsContained[d] + 1 for d in range(dim)]
    # corresponding Points
    extendedPoints = [np.linspace(start[d], end[d], extendedNumPoints[d]) for d in range(dim)]
    interp_mesh = np.array(np.meshgrid(*extendedPoints))
    inter_points = [g.ravel() for g in interp_mesh]  # flattened point mesh
    pointCoordinates = [np.linspace(start[d], end[d], numPoints[d]) for d in range(dim)]
    gridValuesInterpolated = interpn(pointCoordinates, gridValues.reshape(*reversed(numPoints)).transpose(),
                                     list(zip(*inter_points)))
    endTime = time.time()
    # calculate grid combinations
    indices = [list(range(numGrids)) for numGrids in numberOfGridsContained]
    # print indices
    gridCoord = list(zip(*[g.ravel() for g in np.meshgrid(*indices)]))
    # calculate point coordinates
    indices = [list(range(numPerDim)) for numPerDim in numPoints]
    # print indices
    points = list(zip(*[g.ravel() for g in np.meshgrid(*indices)]))
    gridIntegrals = []
    startoffsetGrid = np.zeros(dim)
    endoffsetGrid = np.zeros(dim)
    startTime = time.time()
    for g in range(np.prod(numberOfGridsContained)):
        offsetGrid = np.zeros(dim)
        for d in range(dim):
            offsetGrid[d] = (numPoints[d] - 1) * float(gridCoord[g][d])
            startoffsetGrid[d] = start[d] + float(gridCoord[g][d]) / float(numberOfGridsContained[d]) * length[d]
            endoffsetGrid[d] = start[d] + float(gridCoord[g][d] + 1) / float(numberOfGridsContained[d]) * length[d]
        # copy part of interpolated grid to gridsize
        for i in range(gridsize):
            position = np.inner(extendedOffsets, np.array(points[i]) + offsetGrid)
            gridValues[i] = gridValuesInterpolated[int(position)]
        # calculates integral for current subGrid
        currentSliceSize = gridsize
        for d in reversed(list(range(dim))):
            currentSliceSize = int(currentSliceSize / int(numPoints[d]))
            for i in range(currentSliceSize):
                lineValues = np.zeros(int(numPoints[d]))
                for j in range(int(numPoints[d])):
                    lineValues[j] = gridValues[i + j * offsets[d]]
                gridValues[i] = np.trapz(lineValues, dx=extendedSpacing[d])
        gridIntegrals.append((gridValues[0], np.array(startoffsetGrid), np.array(endoffsetGrid)))
    endTime = time.time()
    return gridIntegrals
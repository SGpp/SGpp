import pysgpp
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib
matplotlib.use("TkAgg")

sys.path.append('/home/rehmemk/git/SGpp/MR_Python/Extended')  # nopep8
import scalarFunctions  # nopep8
from scalarFunctions import objFuncSGpp  # nopep8
sys.path.append('/home/rehmemk/git/SGpp/MR_Python/Vector')  # nopep8
from vectorFunctions import vectorObjFuncSGpp   # nopep8
import vectorFunctions  # nopep8
sys.path.append('/ home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri')  # nopep 8
from sgppOkushiri import maxOkushiri1Out  # nopep8

dim = 1
level = 9
pyFunc = scalarFunctions.getFunction('maxOkushiri', dim)
lb = pysgpp.DataVector([0.5])
ub = pysgpp.DataVector([1.5])
refineType = 'regular'
gridType = 'nakbsplineboundary'
maxLevel = level
degrees = [1, 3, 5]
fig = plt.figure(figsize=plt.figaspect(1./len(degrees)))
for k, degree in enumerate(degrees):
    objFunc = objFuncSGpp(pyFunc)
    reSurf = pysgpp.SparseGridResponseSurfaceBspline(
        objFunc, lb, ub, pysgpp.Grid.stringToGridType(gridType), degree)
    reSurf.regular(maxLevel)
    coeffs = reSurf.getCoefficients()
    grid = reSurf.getGrid()
    gridStorage = grid.getStorage()
    levelWiseAlphaNorms = np.zeros(level+1)
    for i in range(grid.getSize()):
        l = gridStorage.getPointLevel(i, 0)
        x = gridStorage.getPointCoordinate(i, 0)
        # print('{}   {}'.format(l, x))
        levelWiseAlphaNorms[l] += coeffs[i]**2

    plt.subplot(1, 3, k+1)
    plt.plot(range(len(levelWiseAlphaNorms)), levelWiseAlphaNorms)
    plt.gca().set_yscale('log')
    plt.title(degree)
    plt.xlabel('level')
    plt.ylabel('||alpha||')
pyFunc.cleanUp()
plt.show()

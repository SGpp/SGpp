import pysgpp
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib import cm
matplotlib.use("TkAgg")

sys.path.append('/home/rehmemk/git/SGpp/MR_Python/Extended')  # nopep8
import scalarFunctions  # nopep8
from scalarFunctions import objFuncSGpp  # nopep8
sys.path.append('/ home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri')  # nopep 8

dim = 2
level = 10
pyFunc = scalarFunctions.getFunction('maxOkushiri', 2)
lb = pysgpp.DataVector([0.5, 0.5])
ub = pysgpp.DataVector([1.5, 1.5])
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
    levelWiseAlphaNorms = np.zeros((level+1, level + 1))
    summarizedAlphaNorms = np.zeros(level+1)
    for i in range(grid.getSize()):
        lx = gridStorage.getPointLevel(i, 0)
        ly = gridStorage.getPointLevel(i, 1)
        levelWiseAlphaNorms[lx, ly] += coeffs[i]**2
    print('degree {}:'.format(degree))
    print(levelWiseAlphaNorms)
    print('\n')

    coeffMin = sys.float_info.max
    coeffMax = -sys.float_info.max
    for lx in range(level+1):
        for ly in range(level+1-lx):
            if levelWiseAlphaNorms[lx, ly] < coeffMin:
                coeffMin = levelWiseAlphaNorms[lx, ly]
            if levelWiseAlphaNorms[lx, ly] > coeffMax:
                coeffMax = levelWiseAlphaNorms[lx, ly]
    print("min: {}".format(coeffMin))
    print("max: {}".format(coeffMax))

    viridis = cm.get_cmap('viridis')
    normalize = coeffMax - coeffMin
    ax = plt.subplot(2, 3, k+1)
    for lx in range(level+1):
        for ly in range(level+1-lx):
            summarizedAlphaNorms[lx+ly] += levelWiseAlphaNorms[lx, ly]
            #color = viridis(levelWiseAlphaNorms[lx, ly]/normalize)
            color = viridis((np.log(levelWiseAlphaNorms[lx, ly]) - np.log(
                coeffMin)) / (np.log(coeffMax) - np.log(coeffMin)))
            rect = patches.Rectangle(
                (0.1+lx*0.1, 0.1+ly*0.1), 0.05, 0.05, linewidth=1, edgecolor=color, facecolor=color)
            ax.add_patch(rect)
    ax.set_ylim([0., 0.2+level*0.1])
    ax.set_xlim([0., 0.2+level*0.1])
    ax.axis('off')
    plt.title(degree)
    sm = plt.cm.ScalarMappable(cmap=viridis, norm=colors.LogNorm(
        vmin=coeffMin, vmax=coeffMax))
    plt.colorbar(sm)
    ax = plt.subplot(2, 3, 3+k+1)
    plt.plot(range(level+1), summarizedAlphaNorms)
    ax.set_yscale('log')
    plt.ylabel('level')

plt.show()

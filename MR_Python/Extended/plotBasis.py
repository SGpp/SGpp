import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def plotBasis(basis, l, i):
    x = np.linspace(0, 1, 400)
    B = np.zeros(len(x))
    for j in range(len(x)):
        B[j] = basis.eval(l, i, x[j])
    plt.plot(x, B)


def plotWholeLevel(basis, level):
    # # hierarchical
    indices = range(1, (2 ** level - 1) + 1, 2)
    # # all
    # indices = range(1, (2 ** level - 1) + 1)
    print(indices)
    for i in indices:
        plotBasis(basis, level, i)
    points = [i / (2.0 ** level) for i in range(0, 2 ** level + 1)]
    plt.xticks(points)

        
def plotSeveralLevels(basis, levels):
    s = 1
    for l in levels:
        ax = plt.subplot(len(levels), 1, s)
        s += 1
        plotWholeLevel(basis, l)


degree = 3

# plt.figure()
# basis = pysgpp.SNakBsplineBase(degree)
# plotSeveralLevels(basis, [1, 2, 3, 4])
#  
# plt.figure()
# basis = pysgpp.SNakBsplineModifiedBase(degree)
# plotSeveralLevels(basis, [1, 2, 3, 4])
#  
plt.figure()
basis = pysgpp.SNakBsplineExtendedBase(degree)
plotSeveralLevels(basis, [1, 2, 3, 4])

plt.show()

import pysgpp

maxLevel = 10
for dim in range(1, 5):
    print("\ndim: {}".format(dim))
    for l in range(maxLevel+1):
        gridB = pysgpp.Grid_createNakBsplineBoundaryGrid(dim, 3)
        gridB.getGenerator().regular(l)
        grid = pysgpp.Grid_createNakBsplineExtendedGrid(dim, 3)
        grid.getGenerator().regular(l)
        print("level {}     boundary: {} | no boundary: {}".format(l,
                                                                   gridB.getSize(), grid.getSize()))

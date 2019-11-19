import pysgpp
import matplotlib.pyplot as plt

dim = 2
level = 3
gridB = pysgpp.Grid_createNakBsplineBoundaryGrid(dim, 3)
gridB.getGenerator().regular(level)
gridBStorage = gridB.getStorage()
for i in range(gridB.getSize()):
    p = gridBStorage.getPointCoordinates(i)
    plt.plot(p[0],p[1],'b+')
print(gridB.getSize())

grid = pysgpp.Grid_createNakBsplineExtendedGrid(dim, 3)
grid.getGenerator().regular(level)
gridStorage = grid.getStorage()    
plt.figure()
for i in range(grid.getSize()):
    p = gridStorage.getPointCoordinates(i)
    plt.plot(p[0],p[1],'b+')
print(grid.getSize())    

plt.show()

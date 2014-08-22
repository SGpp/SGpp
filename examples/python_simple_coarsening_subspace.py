#!/usr/bin/python

# import modules
import sys
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import *
import matplotlib.pyplot as plotter

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (dim)


# create regular grid, level 3
level = 5
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "Start: number of grid points:  %d" % (gridStorage.size())

# definition of function to interpolate - nonsymmetric(!)
f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1-x1
# create coefficient vector
alpha = DataVector(gridStorage.size())

#store old files
xCoordsOld = []
yCoordsOld = []

for i in xrange(gridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        gridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordsOld.append(gridPointCoordinates[0])
        yCoordsOld.append(gridPointCoordinates[1])

# now refine adaptively 5 times
for refnum in range(5):
    # set function values in alpha
    for i in xrange(gridStorage.size()):
        gp = gridStorage.get(i)
        alpha[i] = f(gp.abs(0), gp.abs(1))
 
    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)
     
     
    #coarsening
    coarsening = SubspaceCoarsening()
    functor = SurplusCoarseningFunctor(alpha,1)
    coarsening.free_coarsen(gridStorage,functor,alpha)
    #decorator.createSubspace(gridStorage,)
   
    print "Refinement step %d, new grid size: %d" % (refnum+1, gridStorage.size())
    
    #
    #plot grid
    #
 
    
    #initialize plotter
    plotter.hold(True)
    
    #plot function
#     X = numpy.arange(0,1,0.1)
#     Y = numpy.arange(0,1,0.1)
#     X,Y = numpy.meshgrid(X, Y)
#     Z = f(X,Y)
#     
#     figure = plotter.fit
#     ax = figure.add_subplot(111,projetion='3d')
#     ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
    
 
    xCoordinates = []
    yCoordinates = []
 
    #print all points
    for i in xrange(gridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        gridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordinates.append(gridPointCoordinates[0])
        yCoordinates.append(gridPointCoordinates[1])
    
    plotter.scatter(xCoordsOld, yCoordsOld, c='r')
    plotter.scatter(xCoordinates, yCoordinates, c='b')
    xCoordsOld = xCoordinates
    yCoordsOld = yCoordinates
 
    #show plot
 
    plotter.hold(False)
    plotter.show()
 
    # extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.size())

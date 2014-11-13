#!/usr/bin/python

# import modules
import sys
from numpy.core.function_base import linspace
from Crypto.Util.number import size
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import *
import matplotlib.pyplot as plotter
import numpy
import csv



def serializeToCSV(path,xvec,yvec):
    
    with open(path, 'w') as csvFile:
            csvWriter = csv.writer(csvFile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            csvWriter.writerow(['x0','x1'])
            
            for i in xrange(1,len(xvec)):
                data = [xvec[i],yvec[i]]  
                csvWriter.writerow(data)

            del csvWriter
            del csvFile



# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createModLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (dim)


# create regular grid, level 3
level = 1
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "Start: number of grid points:  %d" % (gridStorage.size())

# definition of function to interpolate - nonsymmetric(!)
#f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1-x1
f = lambda x0, x1: x0**2+x1**2
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

refinement = HashRefinement()
decorated = RefinementDecorator(refinement)
decorator = SubspaceGSGRefinement(decorated,dim)
# refine a single grid point each time



for refnum in range(15):
    # set function values in alpha
    for i in xrange(gridStorage.size()):
        gp = gridStorage.get(i)
        alpha[i] = f(gp.abs(0), gp.abs(1))
 
    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)
     
    functor = SurplusRefinementFunctor(alpha,1)
    #refinement  stuff
    decorator.freeRefineSubspace(gridStorage,functor)
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
     
    plotter.scatter(xCoordinates, yCoordinates, c='b')
    plotter.scatter(xCoordsOld, yCoordsOld, c='r')
    xCoordsOld = xCoordinates
    yCoordsOld = yCoordinates
 
    #show plot
 
    plotter.hold(False)
    plotter.show()
 
    # extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.size())

serializeToCSV("GSGMean_Sum.csv",xCoordsOld,yCoordsOld)
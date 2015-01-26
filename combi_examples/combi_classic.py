# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

'''
Python example for the combination technique. Be sure that the folder containing
pysgpp.py can be found in the PYTHONPATH
'''

from pysgpp import SerialCombiGrid,S_CT,DoubleVector

def f_3D(coords):
    return 1.0 + (0.25 * (coords[0] - 0.7) * (coords[0] - 0.7) + 2.0)\
            + (0.25 * (coords[1] - 0.7) * (coords[1] - 0.7) + 2.0)\
            + (0.25 * (coords[2] - 0.7) * (coords[2] - 0.7) + 2.0)



def main():
    dim=3
    level=3
    
    #setup combischeme
    scheme=S_CT(dim,level)
    
    #print subspaces
    levels=scheme.getLevels()
    coeffs=scheme.getCoef()
    for level,coef in zip(levels,coeffs):
        print level,coef
    
    #setup grid
    grid=SerialCombiGrid(scheme,True)
    
    #create full grids
    grid.createFullGrids()
    
    #fill data
    for i in range(grid.getNrFullGrid()):
        fgrid=grid.getFullGrid(i)
        data=fgrid.getElementVector()
        coord=DoubleVector(dim)
        for j in range(fgrid.getNrElements()):
            fgrid.getCoords(j,coord)
            data[j]=f_3D(coord)
    
    #evaluate the grid
    eval_coordinates=DoubleVector([0.1,0.1,0.1])
    print 'function value:\t',f_3D(eval_coordinates)
    print 'combined value:\t',grid.eval(eval_coordinates)

if __name__ == '__main__':
    main()
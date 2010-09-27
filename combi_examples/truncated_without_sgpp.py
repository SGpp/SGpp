#!/usr/bin/python
#/*
 #* clasical_combi.cpp
 #*
 #*  Created on: Aug 20, 2010
 #*      Author: Aliz Nagy
 #*
 #*      A program presenting some examples with the basic functions introduced with the combination technique and their usage without actually using a sparse grid 
 #      The combination of the fullgrids in this case adds up to a modified linear trapezoid boundary grid
 #*	The code includes the construction of a fullgridset of given level and given dimension, the assignment of function values to the points of the fullgrid,
 #* 	the evaluation of the fullgrids and the combination of the result, and the comparison with the actual function value
 #*/
# import modules
import sys
import re
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import DataVector, Grid, FullGrid, FullGridSet
from ctypes import *


dim = 3
level = 4
l_user = 2
#The function we want to interpolate
f = lambda x0,x1,x2: 1.0+(0.25*(x0-0.7)*(x0-0.7)+2.0)+(0.25*(x1-0.7)*(x1-0.7)+2.0)+(0.25*(x2-0.7)*(x2-0.7)+2.0)
# Create the set of fullgrids, the gridType can be "linear", "linearBoundary", "linearTrapezoidBoundary or "squareRoot"
fgs = FullGridSet(dim,level,l_user)
print "Number of grids",fgs.getSize()
print "The grids:"
#This prints the levels of all fullgrids
fgs.printGridSet()


#Create a new datavector which contains the coordinates of the point we want to interpolate
p=DataVector(dim);
p[0] = 0.91;
p[1] = 0.23;
p[2] = 0.76;
#Fill the full grids with the function values
for i in xrange(fgs.getSize()):
    fg=fgs.at(i)  
    m=fg.getSize()
    for j in xrange(m):
       #Set the value
        fg.set(j,f(fg.getCoord(0,j),fg.getCoord(1,j),fg.getCoord(2,j)))      
       # Evaluates the fullgrid in an arbitrary point, and assigns the resulting value to the field variable value of the fullgrid,
       # The same value is returned by the function and can be accesed later through the function call fullgrid.val()
        fg.eval(p)
 #Combines the interpolated results on the fullgrids into one value, which in case of function interpolation equals the value on the sparse grid itself
 #We print the interpolation value as Uc(p)
 #This and the evals on fullgrids can be replaced by the eval() function of the fullgrids
print "Uc(p)=",fgs.combinedResult()
#We will now verify if the same value is obtained on the sparse grid

#We print the real value of the function in the given point

print "f(p)=",f(p[0],p[1],p[2])

#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_dataVectorSerializeDemo_py Using the DataVector object
## 
## This example shows how to initialize a DataVector object, store it to a file
## and then to restore it back.

from pysgpp import *

## Create a vector with 3 elements
v = DataVector(3)
v.append(1.0)
v.append(2.0)
v.append(3.0)

## Store the vector to file
print( v )
v.toFile("dataVectorTest.mat")

## Load the created DataVector from file and then save it again (for no particular reason)
v2 = DataVector.fromFile("dataVectorTest.mat")
print( v2 )

v2.toFile("dataVectorTest2.mat")

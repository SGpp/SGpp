#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_dataMatrixSerializeDemo_py Using the DataMatrix object
## 
## This example shows how to initialize a DataMatrix object, store it to a file
## and then to restore it back.

from pysgpp import *

## Create a 2x2 matrix
m = DataMatrix(2,2)
m.set(0, 0, 1.0)
m.set(0, 1, 2.0)
m.set(1, 0, 3.0)
m.set(1, 1, 4.0)

## Store the matrix to file
m.toFile("dataMatrixTest.mat")

## Load the created DataMatrix from file and then save it again (for no particular reason)
m2 = DataMatrix.fromFile("dataMatrixTest.mat")

m2.toFile("dataMatrixTest2.mat")

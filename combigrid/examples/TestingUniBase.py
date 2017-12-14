#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_bSplines_py bSplines.py
## plots anisotropic full grids that form part of the combination technique

from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation,\
    CombigridMultiOperation, CombigridTensorOperation
import pysgpp
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, dblquad
from pysgpp.extensions.datadriven.uq.dists import Uniform
from pysgpp.extensions.datadriven.uq.dists.Beta import Beta
from numpy import square

x = 0.337
degree = 3
level = 1
index = 2

points=[1]*((2**level)+1)
for i in range(0,(2**level)+1):
    points[i] = (i*1.0)/(2**level)

evalpoints = np.linspace(0,1,200)

evNew=[0]*len(evalpoints)
ev=[0]*len(evalpoints)

xi = pysgpp.DoubleVector()
pysgpp.createdeg3NakKnots(points, xi)
for i in range(0,len(evalpoints)):
    evNew[i] = pysgpp.expUniformNaKBspline(evalpoints[i], degree, index, points)
    ev[i]  = pysgpp.LagrangePolynomial(evalpoints[i],points,index)
    #ev[i] = pysgpp.nonUniformBSpline(evalpoints[i],degree,index,xi)
    
print(points)
print(min([x1 - x2 for (x1, x2) in zip(evNew, ev)]))
    
fig, ax = plt.subplots()
plt.plot(evalpoints,evNew,label='New',marker='*')
plt.plot(evalpoints,ev,label='Old')
legend = ax.legend(loc='upper left')
for label in legend.get_texts():
    label.set_fontsize('x-small')
plt.show()


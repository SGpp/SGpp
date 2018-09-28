#!/usr/bin/python
#Copyright(C) 2008 - today The SG++ project
#This file is part of the SG++ project.For conditions of distribution and
#use, please see the copyright notice provided with SG++ or at
#sgpp.sparsegrids.org

## \page example_point_distributions_py Point Distributions (Python)
## This simple example demonstrates the different types of 1-D point distributions
## available in the combigrid module.

## @image html point_distributions.png ""

import numpy as np
import matplotlib.pyplot as plt

import pysgpp

n_points = 13
y_range = np.ones([1, n_points])
points = np.zeros([6, n_points])

uni = pysgpp.UniformPointDistribution()
unib = pysgpp.UniformBoundaryPointDistribution()
cheby = pysgpp.ChebyshevDistribution()
ccurtis = pysgpp.ClenshawCurtisDistribution()
leja = pysgpp.LejaPointDistribution()
l2_leja = pysgpp.L2LejaPointDistribution()


for i in range(n_points):
    points[0, i] = uni.compute(n_points, i)
    points[1, i] = unib.compute(n_points, i)
    points[2, i] = cheby.compute(n_points, i)
    points[3, i] = ccurtis.compute(n_points, i)
    points[4, i] = leja.compute(n_points, i)
    points[5, i] = l2_leja.compute(n_points, i)
    
plt.scatter(points[0,], y_range*6, label="Uniform")
plt.scatter(points[1,], y_range*5, label="UniformBoundary")
plt.scatter(points[2,], y_range*4, label="Chebychev")
plt.scatter(points[3,], y_range*3, label="Clenshaw-Curtis")
plt.scatter(points[4,], y_range*2, label="Leja")
plt.scatter(points[5,], y_range*1, label="L2-Leja")

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
           ncol=3, fancybox=True)

plt.savefig("point_distributions.png")

plt.show()


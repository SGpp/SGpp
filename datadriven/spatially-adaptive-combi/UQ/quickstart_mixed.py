##############################################################################
#
# This is a simple quickstart python program to show how to use chaospy
# to quantify the uncertainty of a simple function with a discontinuity.
#
# The goal is to create a sparse adaptive surrogate that can also be used
# instead of the model() function.
#
# Uncertain parameter:
# - parameter1
# - parameter2
# - parameter3
#
# Quantity of interest:
# - some value
#
# Author: Florian Kuenzner
#
##############################################################################


import chaospy as cp
import numpy as np
import math
import json
import os
from sys import path
path.append('../')
from Function import *
from spatiallyAdaptiveExtendSplit import *
from ErrorCalculator import *

#################################################################################################
# parameter setup
parameter1 = 0.3 # 0.3
parameter1_var = 0.03
parameter1_min = 0.1 #0.1
parameter1_max = 0.3 #0.5
parameter2 = 1.0
parameter2_var = 0.5 #0.03
parameter2_min = 1.0 #0.8
parameter2_max = 1.2 #1.2
parameter3 = 1.6
parameter3_var = 0.3
parameter3_min = 1.4 #1.4
parameter3_max = 1.8 #1.8

#################################################################################################
# setup uncertain parameter
#parameter1Dist = cp.Uniform(parameter1_min, parameter1_max)
parameter1Dist = cp.Normal(parameter1, parameter1_var)
#parameter2Dist = cp.Uniform(parameter2_min, parameter2_max)
parameter2Dist = cp.Normal(parameter2, parameter2_var)
#parameter3Dist = cp.Uniform(parameter3_min, parameter3_max)
parameter3Dist = cp.Normal(parameter3, parameter3_var)
dist = cp.J(parameter1Dist, parameter2Dist, parameter3Dist)

#################################################################################################
# generate nodes and weights
#q = 3  # number of collocation points for each dimension
#nodes, weights = cp.generate_quadrature(q, dist, rule="G")
mean=[parameter1,parameter2,parameter3]
std_dev=[parameter1_var, parameter2_var, parameter3_var]
a = [0.2,-1.3,0]
b = [0.4,3,3]
# weight function that applies the normal distribution in x and y and leaves z unchanged (uniform distributed)
weight_function = lambda x : np.prod([norm.pdf(x=x[d], loc=mean[d], scale=std_dev[d]) / (norm.cdf(b[d],loc=mean[d], scale=std_dev[d]) - norm.cdf(a[d],loc=mean[d], scale=std_dev[d])) for d in range(2)])
model = FunctionUQWeighted(FunctionUQ(), weight_function)
gridX = TruncatedNormalDistributionGrid1D(a=a[0],b=b[0],mean=parameter1,std_dev=parameter1_var)
gridY = TruncatedNormalDistributionGrid1D(a=a[1],b=b[1],mean=parameter2,std_dev=parameter2_var)
gridZ = GaussLegendreGrid1D(a=a[2],b=b[2])
grid = MixedGrid(a=a,b=b,dim=3,grids=[gridX,gridY,gridZ])
errorOperator2=ErrorCalculatorExtendSplit()
adaptiveCombiInstanceExtend = SpatiallyAdaptiveExtendScheme(a, b,2,grid,version=0)
adaptiveCombiInstanceExtend.performSpatiallyAdaptiv(1,2,model,errorOperator2,10**-10, do_plot=False)
nodes, weights = adaptiveCombiInstanceExtend.get_points_and_weights()
print("Number of points:", len(nodes))
print("Sum of weights:", sum(weights))
weights = np.asarray(weights) * 1.0/sum(weights)
nodes_transpose = list(zip(*nodes))

#################################################################################################
# propagate the uncertainty
value_of_interests = [model(node) for node in nodes]
value_of_interests = np.asarray(value_of_interests)
print("Mean", np.inner(weights, value_of_interests))
#################################################################################################
# generate orthogonal polynomials for the distribution
OP = cp.orth_ttr(3, dist)

#################################################################################################
# generate the general polynomial chaos expansion polynomial
gPCE = cp.fit_quadrature(OP, nodes_transpose, weights, value_of_interests)

#################################################################################################
# calculate statistics
E = cp.E(gPCE, dist)
StdDev = cp.Std(gPCE, dist)

#print the stastics
print("mean: %f" % E)
print("stddev: %f" % StdDev)

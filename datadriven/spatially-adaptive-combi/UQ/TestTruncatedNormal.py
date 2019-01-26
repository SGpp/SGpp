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
model = FunctionUQ()
#a = [parameter1_min, parameter2_min, parameter3_min]
#b = [parameter1_max, parameter2_max, parameter3_max]
#grid = TruncatedNormalDistributionGrid(mean=[2], std_dev=[0.5])
#grid.setCurrentArea([0],[100],[3])
grid = TruncatedNormalDistributionGrid(mean=[0.3], std_dev=[0.03], global_a=[-100],global_b=[100])
grid.setCurrentArea([-100],[100],[3])
print("grid",grid.get_points_and_weights())
points, weights = grid.get_points_and_weights()
print(points[0][0])
print(weights)
print(sum([math.sin(float(points[d][0])) * float(weights[d]) for d in range(len(points))]))


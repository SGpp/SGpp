#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_bSplines_py bSplines.py
## plots anisotropic full grids that form part of the combination technique

try:
    from argparse import ArgumentParser
    from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
    from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
    from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
    from pysgpp.pysgpp_swig import DataVector, CombigridOperation,\
        CombigridMultiOperation, CombigridTensorOperation
    from pysgpp.extensions.datadriven.uq.dists import Uniform
    from pysgpp.extensions.datadriven.uq.dists.Beta import Beta

    import pysgpp
    import os

    import matplotlib.pyplot as plt
    import numpy as np

    from numpy import square
    from scipy.integrate import quad, dblquad

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)


# calulate scalar product of basis functions
#===============================================================================
# For basis functions (1,1) (0,0) (0,1) (2,1) (2,3) in this order (the shift between level 0 and 1 is due to the combigrid conversion
#
# 1                  0.6666666666666666  0.5                0.275                0.275
# 0.6666666666666666 0.5333333333333333  0.25               0.2504166666666666   0.11041666666666665
# 0.5                0.25                0.3333333333333333 0.0675               0.2075
# 0.275              0.2504166666666666  0.0675             0.13857142857142857  0.016857142857142862
# 0.275              0.11041666666666665 0.2075             0.016857142857142862 0.13857142857142857
#===============================================================================

degree = 5
basis = pysgpp.SNakBsplineBoundaryCombigridBase(degree)
l1 = 1
i1 = 1
l2 = 1
i2 = 1
scalarProduct = quad(lambda x: basis.eval(l1,i1,x) * basis.eval(l2,i2,x), 0,1)
print(scalarProduct)

# print basis functions
B00=[];B01=[]
B11=[]
B21=[];B23=[]
B31=[];B33=[];B35=[];B37=[]
B41=[];B43=[];B45=[];B47=[];B49=[];B411=[];B413=[];B415=[]
B51=[];B53=[];B55=[];B57=[];B59=[];B511=[];B513=[];B515=[];B517=[];B519=[];B521=[];B523=[];B525=[];B527=[];B529=[];B531=[]
X = np.linspace(0,1,250)
for i in range(0,len(X)):
    B00.append(basis.eval(0,0,X[i]));    B01.append(basis.eval(0,1,X[i]));
    B11.append(basis.eval(1,1,X[i]))
    B21.append(basis.eval(2,1,X[i]));    B23.append(basis.eval(2,3,X[i]))
    B31.append(basis.eval(3,1,X[i]));    B33.append(basis.eval(3,3,X[i]));    B35.append(basis.eval(3,5,X[i]));    B37.append(basis.eval(3,7,X[i]))
    B41.append(basis.eval(4,1,X[i]));    B43.append(basis.eval(4,3,X[i]));    B45.append(basis.eval(4,5,X[i]));    B47.append(basis.eval(4,7,X[i]));B49.append(basis.eval(4,9,X[i]));    B411.append(basis.eval(4,11,X[i]));    B413.append(basis.eval(4,13,X[i]));    B415.append(basis.eval(4,15,X[i]))
    B51.append(basis.eval(5,1,X[i]));    B53.append(basis.eval(5,3,X[i]));    B55.append(basis.eval(5,5,X[i]));    B57.append(basis.eval(5,7,X[i]));B59.append(basis.eval(5,9,X[i]));    B511.append(basis.eval(5,11,X[i]));    B513.append(basis.eval(5,13,X[i]));    B515.append(basis.eval(5,15,X[i])); B517.append(basis.eval(5,17,X[i]));    B519.append(basis.eval(5,19,X[i]));    B521.append(basis.eval(5,21,X[i]));    B523.append(basis.eval(5,23,X[i]));B525.append(basis.eval(5,25,X[i]));    B527.append(basis.eval(5,27,X[i]));    B529.append(basis.eval(5,29,X[i]));    B531.append(basis.eval(5,31,X[i]))
fig = plt.figure()
#plt.plot(X,B00); plt.plot(X,B01)
#plt.plot(X,B11)
#plt.plot(X,B21);plt.plot(X,B23)
#plt.plot(X,B31);plt.plot(X,B33);plt.plot(X,B35);plt.plot(X,B37)
plt.plot(X,B41);plt.plot(X,B43);plt.plot(X,B45);plt.plot(X,B47);plt.plot(X,B49);plt.plot(X,B411);plt.plot(X,B413);plt.plot(X,B415)
#plt.plot(X,B51);plt.plot(X,B53);plt.plot(X,B55);plt.plot(X,B57);plt.plot(X,B59);plt.plot(X,B511);plt.plot(X,B513);plt.plot(X,B515);plt.plot(X,B517);plt.plot(X,B519);plt.plot(X,B521);plt.plot(X,B523);plt.plot(X,B525);plt.plot(X,B527);plt.plot(X,B529);plt.plot(X,B531)
plt.xticks(np.linspace(0, 1, 2**4 + 1))
plt.show()

from __future__ import division
import pysgpp
import numpy as np
import matplotlib.pyplot as plt
import math
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d

pi = math.pi

class interpolation_function(object):
  def __init__(self, d, f):
    self.f = f
    self.d = d
    self.grid = pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3)
    self.gridStorage = self.grid.getStorage()
    try :
      self.hierarch = pysgpp.createOperationHierarchisation(self.grid)
    except :
      self.hierarch = pysgpp.createOperationMultipleHierarchisation(self.grid)
    self.opeval = pysgpp.createOperationEvalNaive(self.grid)
    self.alpha = pysgpp.DataVector(self.gridStorage.getSize())

  def create_interpolation(self, grid_lvl):
    self.gridStorage.clear()
    self.grid.getGenerator().regular(grid_lvl)
    self.alpha = pysgpp.DataVector(self.gridStorage.getSize())
    self.min_f = float('inf')
    for i in range(self.gridStorage.getSize()):
      gp = self.gridStorage.getPoint(i)
      x = [self.gridStorage.getCoordinate(gp, j) for j in range(self.d)]
      self.alpha[i] = self.f(x)
      if self.alpha[i] < self.min_f:
          self.min_f = self.alpha[i]
          self.min_x = x
    self.hierarch.doHierarchisation(self.alpha)

  def __call__(self, x):
    if (self.d == 1 and not isinstance(x, list)):
        x = [x]
    return self.opeval.eval(self.alpha, pysgpp.DataVector(x))

# ---------------------------------------------

def distrib(x):
    if isinstance(x, list):
        x = x[0]
    return (np.sin(2.5*pi*x - (pi - pi / 4.)) + 1. / 2**0.5) / 0.527044

def parabola(x):
  res = 1.
  for i in range(len(x)):
    res *= x[i] * (1. - x[i]) * 4.;
  return res

def const_0(x):
  return 0.0

def eval_rosenblatt1d(sg_pdf, xs):
    op = pysgpp.createOperationRosenblattTransformation1D(sg_pdf.grid)
    ys = []
    for i,x in enumerate(xs):
        print("---------------{}---------------".format(i))
        ys.append(op.doTransformation1D(sg_pdf.alpha, x))
    return ys
  # ---------------------------------------------

def eval_rosenblattdd(sg_pdf, xs):
    op = pysgpp.createOperationRosenblattTransformation(sg_pdf.grid)
    X, Y = np.meshgrid(xs, xs)
    input_points = pysgpp.DataMatrix(len(xs), sg_pdf.d)
    output_points = pysgpp.DataMatrix(len(xs), sg_pdf.d)
    for i in range(len(xs)):
      for j in range(sg_pdf.d):
        input_points.set(i, j , 0.5)
    op.doTransformation(sg_pdf.alpha, input_points, output_points)
    print(output_points)

  # ---------------------------------------------


def eval_inverse_rosenblatt1d(sg_pdf, xs):
  op = pysgpp.createOperationInverseRosenblattTransformation1D(sg_pdf.grid)
  ys = []
  print(op.doTransformation1D(sg_pdf.alpha, 0.999942))
  for i,x in enumerate(xs):
    # print("---------------{}---------------".format(i))
    ys.append(op.doTransformation1D(sg_pdf.alpha, x))
  return ys

xs = np.arange(0., 1.01, 0.01)
l_max = 2
d = 2
interpolation = interpolation_function(d, parabola)
interpolation.create_interpolation(l_max)

# test()
# grid_points = np.arange(0, 1.01, 2**-l_max)
# ys = [interpolation(x) for x in xs]

# plotSG2d(interpolation.grid, interpolation.alpha)
# ys = eval_inverse_rosenblatt1d(interpolation, xs)
# ys = eval_rosenblatt1d(interpolation, xs)
# print(ys)
eval_rosenblattdd(interpolation, xs)
# plt.plot(xs, ys)
# plt.scatter(grid_points, np.zeros_like(grid_points))
# plt.pcolormesh(X, Y, Z)
# plt.show()

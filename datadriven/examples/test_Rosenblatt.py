import pysgpp
import numpy as np
import matplotlib.pyplot as plt
import math

pi = math.pi

class interpolation_function():
  def __init__(self, d, f):
    self.f = f
    self.d = d
    self.grid = pysgpp.Grid.createPolyGrid(d, 3)
    self.gridStorage = self.grid.getStorage()
    self.hierarch = pysgpp.createOperationHierarchisation(self.grid)
    self.opeval = pysgpp.createOperationEval(self.grid)
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
    print(self.alpha)
    self.hierarch.doHierarchisation(self.alpha)
    print(self.alpha)

  def __call__(self, x):
    if (self.d == 1 and not isinstance(x, list)):
        x = [x]
    return self.opeval.eval(self.alpha, pysgpp.DataVector(x))

# ---------------------------------------------

def distrib(x):
    if isinstance(x, list):
        x = x[0]
    return (np.sin(2.5*pi*x - (pi - pi/4.)) + 1./2**0.5) / 0.527044

def eval_rosenblatt(sg_pdf, xs):
    op = pysgpp.createOperationRosenblattTransformation1D(sg_pdf.grid)
    ys = []
    for x in xs:
        ys.append(op.doTransformation1D(sg_pdf.alpha, x))
        print("------------------------------")
    return ys
  # ---------------------------------------------

def eval_inverse_rosenblatt(sg_pdf, xs):
  op = pysgpp.createOperationInverseRosenblattTransformation1D(sg_pdf.grid)
  ys = []
  for x in xs:
    ys.append(op.doTransformation1D(sg_pdf.alpha, x))
    print("------------------------------")
  return ys

l_max = 4
interpolation = interpolation_function(1, distrib)
interpolation.create_interpolation(l_max)

xs = np.arange(0, 1.01, 0.01)
grid_points = np.arange(0, 1.01, 2**-l_max)
# ys = [interpolation(x) for x in xs]
ys = eval_rosenblatt(interpolation, xs)

grid_points = np.arange(0, 1.01, 2**-l_max)
grid_point_values = eval_rosenblatt(interpolation, grid_points)

diffs = []


def test():
  ys = [0, 0.153888, 0.307776, 0.461665, 0.615553, 1.52762, 2.39577, 3.01498, 3.23902, 3.01498, 2.39577, 1.52762, 0.615553, 0.461665, 0.307776, 0.153888, 0]

  xs = [i*2**-4 for i in range(0, 2**4 + 1)]
  plt.plot(xs ,ys)

for i in range(1, len(grid_points)):
  diffs.append(grid_point_values[i] - grid_point_values[i-1])
# print(xs)
# print(ys)
for i in range(1, len(ys)):
  if (ys[i] <= ys[i-1]):
    print("ERROR {} <= {}".format(ys[i], ys[i-1]))
print(grid_points)
print(diffs)
plt.plot(xs, ys)
# test()
plt.scatter(grid_points, np.zeros_like(grid_points))
# plt.legend()
plt.show()

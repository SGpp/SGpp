from __future__ import division
from __future__ import print_function
import pysgpp
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis
import numpy as np
import random
import matplotlib.pyplot as plt

def test_base():
    n = 100000
    b = pysgpp.SLinearModifiedBase()
    # print(sum([b.eval(0,0,x) for x in np.linspace(0, 1, n)]) / n)
    for l in range(1,3,1):
        print(l)
        for i in range(1, 2**l, 1):
            s = 1
            s = sum([b.eval(l,i,x) for x in np.linspace(0, 1, n)]) / n
            print(("Level:",l , "index:", i))
            epsilon = abs(s - b.getIntegral(l,i))
            if(epsilon > 10**-6):
                print("s: {} res: {}".format(s, b.getIntegral(l,i)))
                print("error: {}".format(epsilon))

def test_LTwoDot(grid, l):
    res = 10000
    b = grid.getBasis();
    grid.getGenerator().regular(l)
    gridStorage = grid.getStorage()
    size = gridStorage.getSize()
    # print(size)
    m = pysgpp.DataMatrix(size, size)
    opMatrix = pysgpp.createOperationLTwoDotExplicit(m, grid)
    # print m
    # m_comp = pysgpp.DataMatrix(size, size)
    for i in range(gridStorage.getSize()):
        for j in range(i, gridStorage.getSize()):
            gpi = gridStorage.getPoint(i)
            gpj = gridStorage.getPoint(j)
            sol = 1
            # print "--------"
            # print "i:{} j:{}".format(i, j)
            for k in range(d):
                lik = gpi.getLevel(k)
                iik = gpi.getIndex(k)
                ljk = gpj.getLevel(k)
                ijk = gpj.getIndex(k)
                # print "i l,i: {},{}   j l,i: {},{}".format(lik, iik, ljk, ijk)
                xs = np.linspace(0, 1, res)
                tmp = sum([b.eval(lik, iik, x)*b.eval(ljk, ijk, x) for x in xs]) / res
                sol *= tmp
                # print("lik:{} iik:{} ljk:{} ijk:{} k:{} tmp: {}".format(lik, iik, ljk, ijk, k,tmp))
            # print(sol)
            error = abs(m.get(i,j) - sol)
            # print error
            if(error >= 10**-4):
                print("i:{} j:{} error: {}".format(i, j, error))
                print("iik:{} lik:{} ijk:{} ljk:{} error: {}".format(iik, lik, ijk, ljk, error))
                print("is:{} should:{}".format(m.get(i,j), sol))

def test_LTwoDotImplicit(grid, l):
    grid.getGenerator().regular(l)
    gridStorage = grid.getStorage()
    size = gridStorage.getSize()
    # print(size)
    m = pysgpp.DataMatrix(size, size)
    opExplicit = pysgpp.createOperationLTwoDotExplicit(m, grid)
    op = pysgpp.createOperationLTwoDotProduct(grid)
    alpha = pysgpp.DataVector(size)
    resultExplicit = pysgpp.DataVector(size)
    result = pysgpp.DataVector(size)
    for i in range(size):
        alpha[i] = 1
    opExplicit.mult(alpha, resultExplicit)
    op.mult(alpha, result)
    for i in range(size):
        if result[i] != resultExplicit[i]:
            print("Error result entry {} differs".format(i))

        if abs(result[i] - resultExplicit[i]) > 1e-16:
            # print result[i] - resultExplicit[i]
            print("result:{}".format(result[i]))
            print("resultExplicit:{}".format(resultExplicit[i]))

def test_laplace(grid, lmax):
    resolution = 100000
    grid.getGenerator().regular(lmax)
    gridStorage = grid.getStorage()
    size = gridStorage.getSize()
    b = getBasis(grid)
    op = pysgpp.createOperationLaplace(grid)
    alpha = pysgpp.DataVector(size)
    result = pysgpp.DataVector(size)

    for point_i in range(size):
      for point_j in range(size):
        gp_i = gridStorage.getPoint(point_i)
        gp_j = gridStorage.getPoint(point_j)
        print("--------")
        for i in range(0, size):
          alpha[i] = 0
        alpha[point_i] = 1
        op.mult(alpha, result)
        xs = np.linspace(0, 1, resolution)
        approx = sum([b.evalDx(gp_i.getLevel(0), gp_i.getIndex(0), x) * b.evalDx(gp_j.getLevel(0), gp_j.getIndex(0), x) for x in xs]) / resolution
        print("i,j: {},{} result: {} approx:{}".format(point_i, point_j, result[point_j], approx))
        if(abs(result.get(point_j) - approx) > 1e-1):
          print("--------")
          print("points: {},{} ".format(point_i, point_j))
          print("approx:{}".format(approx))
          print("result:{}".format(result.get(point_j)))
          # print result
          print("--------")

def test_poly_evaldx():
    l = 3
    i = 1
    x = 0.12
    eps = 0.0001
    b = pysgpp.SPolyModifiedClenshawCurtisBase(3)
    tang = b.evalDx(l, i, x)
    sec = (b.eval(l, i, x + eps) -  b.eval(l, i, x - eps)) / (2*eps)
    print("evalDx:{}".format(tang))
    print("sekante:{}".format(sec))
    print("evals: {} {}".format( b.eval(l, i, x - eps), b.eval(l, i, x + eps) ))

def plot_evaldx():
    l = 3
    i = 1
    xs = np.linspace(0, 1, 50)
    b = pysgpp.SPolyClenshawCurtisBase(3)
    plt.plot(xs, [b.evalDx(l, i, x) for x in xs])
    plt.show()

def plot_evaldx_prod(grid, lmax, i, j):
    grid.getGenerator().regular(lmax)
    gridStorage = grid.getStorage()
    b = getBasis(grid)
    gp_i = gridStorage.getPoint(i)
    gp_j = gridStorage.getPoint(j)
    lik = gp_i.getLevel(0)
    ljk = gp_j.getLevel(0)
    iik = gp_i.getIndex(0)
    ijk = gp_j.getIndex(0)
    xs = np.linspace(0, 1, 100000)

    plt.plot(xs, [b.evalDx(lik, iik, x) * b.evalDx(ljk, ijk, x) for x in xs])
    plt.show()


def test_laplace2(grid, lmax):
  # grid.getGenerator().regular(lmax)
  gridStorage = grid.getStorage()
  size = gridStorage.getSize()
  m = pysgpp.DataMatrix(size, size)
  op = pysgpp.createOperationLaplaceExplicit(m, grid)
  print(m)



def test_firstMoment(grid, lmax):
  grid.getGenerator().regular(lmax)
  resolution = 100000
  gridStorage = grid.getStorage()
  b = grid.getBasis()
  op = pysgpp.createOperationFirstMoment(grid)
  alpha = pysgpp.DataVector(grid.getSize(), 1.0)
  bounds = pysgpp.DataMatrix(1, 2, 0.0)
  bounds.set(0, 1, 1.0)
  res = 0.0
  for i in range(grid.getSize()):
    lev = gridStorage.getPoint(i).getLevel(0)
    ind = gridStorage.getPoint(i).getIndex(0)
    temp_res = 0.0
    for c in range(resolution):
        x = float(c) / resolution
        temp_res += x * b.eval(lev, ind, x)
    res += alpha.get(i) * temp_res / resolution
  print("--FirstMoment--")
  print(res)
  print(op.doQuadrature(alpha, bounds))
  print (res - op.doQuadrature(alpha, bounds))

def test_secondtMoment(grid, lmax):
  gridStorage = grid.getStorage()
  gridStorage.clear()
  grid.getGenerator().regular(lmax)
  resolution = 1000000
  b = grid.getBasis()
  op = pysgpp.createOperationSecondMoment(grid)
  alpha = pysgpp.DataVector(grid.getSize(), 1.0)
  bounds = pysgpp.DataMatrix(1, 2, 0.0)
  bounds.set(0, 1, 1.0)
  res = 0.0
  for i in range(grid.getSize()):
    lev = gridStorage.getPoint(i).getLevel(0)
    ind = gridStorage.getPoint(i).getIndex(0)
    temp_res = 0.0
    for c in range(resolution):
        x = float(c) / resolution
        temp_res += x * x *b.eval(lev, ind, x)
    res += alpha.get(i) * temp_res / resolution

  print("--SecondMoment--")
  print(res)
  print(op.doQuadrature(alpha, bounds))
  print (res - op.doQuadrature(alpha, bounds))

# test_poly_evaldx()
# plot_evaldx()
# test_base()
d = 1
l = 3
grid = pysgpp.Grid.createModLinearGrid(d)
# plot_evaldx_prod(grid, 4, 1, 4)
# test_laplace(grid, l)
# test_laplace2(grid, l)
# test_LTwoDot(grid, l)
# test_LTwoDotImplicit(grid, l)
test_firstMoment(grid, l)
test_secondtMoment(grid, l)

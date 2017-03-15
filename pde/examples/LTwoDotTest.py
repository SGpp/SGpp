import pysgpp
import numpy as np
import random
import matplotlib.pyplot as plt

def test_base():
    n = 100000
    b = pysgpp.SLinearModifiedBase()
    # print(sum([b.eval(0,0,x) for x in np.linspace(0, 1, n)]) / n)
    for l in range(1,3,1):
        print l
        for i in range(1, 2**l, 1):
            s = 1
            s = sum([b.eval(l,i,x) for x in np.linspace(0, 1, n)]) / n
            print("Level:",l , "index:", i)
            epsilon = abs(s - b.getIntegral(l,i))
            if(epsilon > 10**-6):
                print "s: {} res: {}".format(s, b.getIntegral(l,i))
                print "error: {}".format(epsilon)

def test_LTwoDot(grid, l):
    res = 100000
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
                xs = np.linspace(0,1,res)
                tmp = sum([b.eval(lik, iik, x)*b.eval(ljk, ijk, x) for x in xs])/res
                sol *= tmp
                # print("lik:{} iik:{} ljk:{} ijk:{} k:{} tmp: {}".format(lik, iik, ljk, ijk, k,tmp))
            # print(sol)
            error = abs(m.get(i,j) - sol)
            # print error
            if(error >= 10**-5):
                print "i:{} j:{} error: {}".format(i, j, error)
                print "iik:{} lik:{} ijk:{} ljk:{} error: {}".format(iik, lik, ijk, ljk, error)
                print "is:{} should:{}".format(m.get(i,j), sol)


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
            print "Error result entry {} differs".format(i)

        if abs(result[i] - resultExplicit[i]) > 1e-16:
            # print result[i] - resultExplicit[i]
            print "result:{}".format(result[i])
            print "resultExplicit:{}".format(resultExplicit[i])

def test_poly_evaldx():
    l = 3
    i = 1
    x = 0.12
    eps = 0.0001
    b = pysgpp.SPolyBase(3)
    tang = b.evalDx(l, i, x)
    sec = (b.eval(l, i, x + eps) -  b.eval(l, i, x - eps)) / (2*eps)
    print "evalDx:{}".format(tang)
    print "sekante:{}".format(sec)
    print "evals: {} {}".format( b.eval(l, i, x - eps), b.eval(l, i, x + eps) )
test_poly_evaldx()
# test_base()
# d = 1
# l = 7
# grid = pysgpp.Grid.createBsplineClenshawCurtisGrid(d, 3)
# test_LTwoDot(grid, l)
# test_LTwoDotImplicit(grid, l)

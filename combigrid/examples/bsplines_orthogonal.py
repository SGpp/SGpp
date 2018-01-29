#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# \page example_bSplines_py bSplines.py
# plots anisotropic full grids that form part of the combination technique

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


def f(x):
    return np.exp(-x) * np.cos(4 * x * (1 - x))


def mean():
    return 2. / 3.


def var():
    return 16. / 30. - 4. / 9.


class Basis(object):

    def mean(self):
        return quad(lambda x: self(np.array([[x]])), 0, 1)[0]

    def variance(self):
        return np.var([self(np.array([[x]]))[0] for x in np.random.rand(1000)])
    #         return quad(lambda x: self(np.array([[x]])) ** 2, 0, 1)[0] - e ** 2


class BSplineBasis(Basis):

    def __init__(self, level_index_set, degree):
        self.basis = pysgpp.SNakBsplineBoundaryCombigridBase(degree)
        self.level_index_set = level_index_set
        self.n = len(level_index_set)

    def evaluate(self, level, index, x):
        if level == 1 and index == 1:
            return self.basis.eval(0, 0, x)
        elif level == 0 and index == 0:
            return self.basis.eval(1, 1, x)
        else:
            return self.basis.eval(level, index, x)

    def V(self, xs):
        V = np.ndarray((len(xs), self.n))
        for i, x in enumerate(xs):
            for j, (lj, ij) in enumerate(level_index_set):
                V[i, j] = basis.evaluate(lj, ij, x[0])
        return V

    def G(self):
        G = np.ndarray((self.n, self.n))
        for i, (li, ii) in enumerate(self.level_index_set):
            for j, (lj, ij) in enumerate(self.level_index_set):
                if i <= j:
                    G[i, j] = G[j, i] = quad(lambda x: basis.evaluate(li, ii, x) * basis.evaluate(lj, ij, x), 0, 1)[0]
        return G

    def size(self):
        return self.n


class PolynomialBasis(object):

    def __init__(self, maxDegree):
        config = pysgpp.OrthogonalPolynomialBasis1DConfiguration()
#         config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_LEGENDRE
#         self.basis = pysgpp.OrthogonalPolynomialBasis1D(config)

        config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_JACOBI
        config.polyParameters.alpha_ = 3
        config.polyParameters.beta_ = 5
        self.basis = pysgpp.OrthogonalPolynomialBasis1D(config)

#         self.basis = pysgpp.MonomialFunctionBasis1D()
        self.maxDegree = maxDegree + 1

    def evaluate(self, deg, x):
        return self.basis.evaluate(deg, x)

    def V(self, xs):
        V = np.ndarray((len(xs), self.maxDegree))
        for i, x in enumerate(xs):
            for j in xrange(self.maxDegree):
                V[i, j] = basis.evaluate(j, x[0])
        return V

    def G(self):
        G = np.ndarray((self.maxDegree, self.maxDegree))
        for i in xrange(self.maxDegree):
            for j in xrange(self.maxDegree):
                if i <= j:
                    G[i, j] = G[j, i] = quad(lambda x: basis.evaluate(i, x) * basis.evaluate(j, x), 0, 1)[0]
        return G

    def size(self):
        return self.maxDegree


class LinearCombination(Basis):

    def __init__(self, coeffs, basis, R=None):
        self.coeffs = coeffs
        self.basis = basis
        self.n = self.coeffs.shape[0]

        self.R = R
        if self.R is None:
            self.R = np.identity(self.n)

    def __call__(self, xs):
        # compute Vandermonde matrix
        V = basis.V(xs)

        # apply rotation matrix
        Vr = self.R.dot(V.T).T

        return Vr.dot(self.coeffs)


# -------------------------------------------------------------------------
degree = 3
level = 3

level_index_set = [(0, 0), (0, 1), (1, 1)]
for li in xrange(2, level + 1):
    level_index_set += [(li, ii) for ii in xrange(1, 2 ** li, 2)]

n = len(level_index_set)

# basis = PolynomialBasis(5)
basis = BSplineBasis(level_index_set, degree)

# interpolation matrix and right hand side
n = basis.size()

h = 1. / (n + 1)
grid_points = np.array([np.linspace(h, 1 - h, n, endpoint=True)]).T
rhs = np.array([f(grid_points[i]) for i in xrange(n)])

V = basis.V(grid_points)
coeffs = np.linalg.solve(V, rhs)

surrogate = LinearCombination(coeffs, basis)

# compute the scalar products
G = basis.G()
L = np.linalg.cholesky(G)
R = np.linalg.inv(L)

print "|G - LL^T| = %g" % np.linalg.norm(G - L.dot(L.T))
print "|I - R G R.T| = %g" % np.linalg.norm(np.identity(n) - R.dot(G.dot(R.T)))

Vr = R.dot(V.T).T
ocoeffs = np.linalg.solve(Vr, rhs)

print "|f - (R V^T)^T x| = %g" % np.linalg.norm(rhs - Vr.dot(ocoeffs))

orthogonal_surrogate = LinearCombination(ocoeffs, basis, R=R)

print ocoeffs

print "mean =", mean(), ocoeffs[0][0], surrogate.mean()
print "var =", var(), np.sum(ocoeffs[1:] ** 2), surrogate.variance()

xs = np.array([np.linspace(0, 1, 250)]).T

fig = plt.figure()
plt.plot(xs, [f(xi) for xi in xs], label="f")
plt.plot(xs, surrogate(xs), label="surrogate")
plt.plot(xs, orthogonal_surrogate(xs), label="orthogonal")
plt.scatter(grid_points, np.zeros(grid_points.shape))
plt.vlines(grid_points, np.zeros(grid_points.shape), rhs)
plt.legend()
plt.show()

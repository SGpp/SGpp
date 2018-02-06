#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# \page example_bSplines_py bSplines.py
# plots anisotropic full grids that form part of the combination technique

from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend,\
    load_font_properties
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


def u(x):
    return np.log(np.exp(-x) * np.cos(4 * x * (1 - x)))


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

    def __init__(self, grid_points, degree):
        self.dtype = "bsplines"

        self.degree = degree
        self.evaluator = pysgpp.BSplineInterpolationEvaluator(self.degree)
        self.quadrature = pysgpp.BSplineQuadratureEvaluator(self.degree)
        self.scalarProduct = pysgpp.BSplineScalarProductEvaluator(self.degree)
        self.grid_points = grid_points

        grid_points_double_vector = pysgpp.DoubleVector(grid_points.T[0])
        nak_points_vector = pysgpp.DoubleVector()
        nak_points = []
        for xi in pysgpp.createNakKnots(grid_points_double_vector, self.degree):
            nak_points_vector.push_back(xi)
            nak_points.append(xi)

        self.nak_points = np.array(nak_points)
        self.evaluator.setGridPoints(grid_points_double_vector)
        self.quadrature.setGridPoints(grid_points_double_vector)
        self.scalarProduct.setGridPoints(grid_points_double_vector)

    def V(self, xs):
        V = np.ndarray((len(xs), self.size()))
        for i, x in enumerate(xs):
            self.evaluator.setParameter(pysgpp.FloatScalarVector(x[0]))
            ans = self.evaluator.getBasisValues()  # FloatScalarVectorVector
            for j in xrange(ans.size()):
                V[i, j] = ans[j].getValue()
        return V

    def G(self):
        G = np.ndarray((self.size(), self.size()))
        ans = self.scalarProduct.getBasisValues()

        for i in xrange(ans.size()):
            result = ans[i].getValues()  # FloatScalarVectorVector
            for j in xrange(result.size()):
                G[i, j] = result[j].getValue()
        return G

    def size(self):
        return self.grid_points.shape[0]


class PolynomialBasis(object):

    def __init__(self, maxDegree, dtype="jacobi"):
        config = pysgpp.OrthogonalPolynomialBasis1DConfiguration()

        if dtype == "legendre":
            config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_LEGENDRE
        elif dtype == "jacobi":
            config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_JACOBI
            config.polyParameters.alpha_ = 3
            config.polyParameters.beta_ = 5
        else:
            self.basis = pysgpp.OrthogonalPolynomialBasis1D(config)

        self.basis = pysgpp.OrthogonalPolynomialBasis1D(config)
        self.dtype = dtype

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
            for j in xrange(i, self.maxDegree):
                G[i, j] = G[j, i] = quad(lambda x: basis.evaluate(i, x) * basis.evaluate(j, x), 0, 1)[0]
        return G

    def size(self):
        return self.maxDegree


class LinearCombination(Basis):

    def __init__(self, coeffs, basis, C=None):
        self.coeffs = coeffs
        self.basis = basis
        self.n = self.coeffs.shape[0]

        self.C = C
        if self.C is None:
            self.C = np.identity(self.n)

    def __call__(self, xs):
        # compute Vandermonde matrix
        V = basis.V(xs)

        # apply rotation matrix
        Vr = np.dot(V, self.C.T)

        return Vr.dot(self.coeffs)


# -------------------------------------------------------------------------
degree = 3
level = 5

n = 2**level + 1
grid_points = np.array([np.linspace(0, 1, n, endpoint=True)]).T

for basis in [BSplineBasis(grid_points, degree)]:
    #               PolynomialBasis(5, "jacobi"),
    #               PolynomialBasis(5, "monomial"),
    #               PolynomialBasis(5, "legendre")]:

    # interpolation matrix and right hand side
    V = basis.V(grid_points)
    rhs = np.array([u(grid_points[i][0]) for i in xrange(n)])
    coeffs = np.linalg.solve(V, rhs)

    print "|u - V w| = %g" % np.linalg.norm(rhs - V.dot(coeffs))

    surrogate = LinearCombination(coeffs, basis)

    # compute the scalar products
    G = basis.G()
    L = np.linalg.cholesky(G)
    C = np.linalg.inv(L)

    print "cond(G) =", np.linalg.cond(G)
    print "|C - ((L^T)^-1)^T| =", np.linalg.norm(C - np.linalg.inv(L.T).T)
    print "|G - LL^T| = %g" % np.linalg.norm(G - L.dot(L.T))
    print "|I - C G C^T| = %g" % np.linalg.norm(np.identity(basis.size()) - np.dot(C, np.dot(G, C.T)))

    Vr = np.dot(V, C.T)
    ocoeffs = np.linalg.solve(Vr, rhs)

    print "|u - V C^T w| = %g" % np.linalg.norm(rhs - Vr.dot(ocoeffs))

    orthogonal_surrogate = LinearCombination(ocoeffs, basis, C=C)

    xs = np.array([np.linspace(0, 1, 250)]).T
    u_sg = surrogate(xs)
    u_ortho = orthogonal_surrogate(xs)
    print "|u_sg - u_orth| =", np.sqrt(np.mean((u_sg - u_ortho) ** 2))

    if basis.dtype == "bsplines":
        result = basis.quadrature.getBasisValues()
        mean_ct = 0.0
        for i in xrange(result.size()):
            mean_ct += coeffs[i] * result[i].getValue()

        var_ct = np.sum(ocoeffs ** 2) - mean_ct ** 2
    else:
        mean_ct = ocoeffs[0][0]
        var_ct = np.sum(ocoeffs[1:] ** 2)

    print "-" * 80
    print "mean =", mean_ct, surrogate.mean()
    print "var  =", var_ct, surrogate.variance()

#     fig = plt.figure()
#     plt.plot(xs, [u(xi) for xi in xs], label="u")
#     plt.plot(xs, u_sg, label="surrogate")
#     plt.plot(xs, u_ortho, label="orthogonal")
#     plt.scatter(basis.grid_points, np.zeros(grid_points.shape))
#     plt.vlines(basis.grid_points, np.zeros(grid_points.shape), rhs)
#     plt.title(basis.dtype,
#               fontproperties=load_font_properties())
#     plt.legend()
#
#     plt.show()

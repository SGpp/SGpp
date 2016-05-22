#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
#
"""
@file    hdmr.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 23 12:23:32 2013
@brief ANOVA-(H)igh (D)imensional (M)odel (R)epresentation
decomposition for a Sparse-Grid function

@version  0.1

"""
from pysgpp.extensions.datadriven.uq.estimators.IntegralStrategy import IntegralStrategy
from pysgpp.extensions.datadriven.uq.estimators.MarginalIntegralStrategy import MarginalIntegralStrategy

from pysgpp.extensions.datadriven.uq.operations import (evalSGFunction,
                               isNumerical,
                               discretize,
                               hierarchize,
                               checkInterpolation)
from pysgpp import DataVector

import itertools as it
import numpy as np


class HDMR(object):
    """
    The HDMR class
    """

    def __init__(self, grid, alpha, params, nk=None):
        """
        Constructor
        @param grid: Grid
        @param alpha: surplus vector
        @param params: parameter set
        @param E: expectation value of grid, alpha
        @param V: variance of grid, alpha
        @param nk: maximum length of interactions
        """
        self.__grid = grid
        self.__alpha = alpha
        self.__params = params

        self.__ap = self.__params.activeParams()
        self.__U = self.__ap.getIndependentJointDistribution()
        self.__T = self.__ap.getJointTransformation()
        self.__xlim = self.__U.getBounds()
        self.__dim = len(self.__ap)

        if self.__dim < 2:
            raise AttributeError('dimensionality has to be > 1')

        # check if highest order term is required
        self.__has_highest_order_term = not nk or nk > self.__dim - 1

        if self.__has_highest_order_term:
            self.__nk = self.__dim - 1
        else:
            self.__nk = nk

        self.__expectation_funcs = None
        self.__anova_components = None
        self.__variance_components = None

        self._verbose = True
        self.__marginalization = MarginalIntegralStrategy()
        self.__estimation = IntegralStrategy()
        # self.__estimation = MonteCarloStrategy(npaths=40, isPositive=True)

    def __discretize(self):
        """
        Discretize squared f to obtain a well suited grid for all
        further computations.
        """
        def f(p, val):
            q = self.__T.unitToProbabilistic(p)
            return val ** 2 * self.__U.pdf(q)
        grid, _, _ = discretize(self.__grid, self.__alpha, f,
                                refnums=4, deg=1, epsilon=1e-6)

        # hierarchize without pdf
        gs = grid.getStorage()
        nodalValues = DataVector(gs.size())
        p = DataVector(gs.getDimension())
        for i in xrange(gs.size()):
            gs.get(i).getStandardCoordinates(p)
            nodalValues[i] = evalSGFunction(self.__grid, self.__alpha, p)

        self.__alpha = hierarchize(grid, nodalValues)
        self.__grid = grid

        err = checkInterpolation(self.__grid, self.__alpha, nodalValues)
        if err is True:
            import pdb; pdb.set_trace()

    def __computeMean(self):
        def f(p, val):
            return val
        print "estimate mean:",
        self.__E, _ = IntegralStrategy().estimate(self.__T.vol(),
                                                  self.__grid, self.__alpha,
                                                  f, self.__U, self.__T)

    def __computeVariance(self):
        def f(p, val):
            return (val - self.__E) ** 2
        print "estimate variance:",
        self.__V, _ = IntegralStrategy().estimate(self.__T.vol(),
                                                  self.__grid, self.__alpha,
                                                  f, self.__U, self.__T)

    def getSortedPermutations(self, keys):
        """
        Sort keys with respect (1) to their length and (2) their lexical order
        @param keys:
        """
        ans = [tuple()] * len(keys)
        ix = 0
        for x in sorted(np.unique(map(len, keys))):
            for ck in sorted(filter(lambda k: len(k) == x, keys)):
                ans[ix] = ck
                ix += 1

        return ans

    def __expectation_functions(self):
        """
        Compute the marginalized expectation functions for the ANOVA
        decomposition.
        """
        if self.__nk < 1 or self.__nk > self.__dim - 1:
            raise AttributeError('The truncated order has to be in \
                                  {1, ..., %i}' % (self.__dim - 1,))

        # init
        expec = {}
        U = self.__ap.getIndependentJointDistribution()
        T = self.__ap.getJointTransformation()
        vol = T.vol()

        if self._verbose:
            print "-" * 60

        # add higher order terms
        for k in xrange(self.__nk):
            perms = it.combinations(range(self.__dim), r=k + 1)
            for perm in perms:
                # select dimensions to be integrated
                dd = [d for d in xrange(self.__dim) if d not in perm]

                if self._verbose:
                    print "Explore %s, Integrate: %s" % (perm, dd),

                # -----------------------------------------------
                # Make sure that perm and dd are disjoint sets
                # covering the whole set
                assert sorted(list(perm) + dd) == range(self.__dim)
                assert len(dd) == self.__dim - len(perm)
                assert len(dd) > 0
                # -----------------------------------------------
                # compute first moment

                def f(_, val):
                    return val
                grid, alpha, err = self.__marginalization\
                                       .estimate(vol, self.__grid,
                                                 self.__alpha,
                                                 f, U, T, dd)
                expec[perm] = grid, alpha

                if self._verbose:
                    print "L2 err = %g" % err

        # add highest order term
        if self.__has_highest_order_term:
            perm = tuple(range(self.__dim))
            expec[perm] = self.__grid, self.__alpha

        if self._verbose:
            print "-" * 60

        return expec

    def __combine_terms(self, perm):
        """
        Combine the terms in order to get the ANOVA component specified in perm
        @param perm: tuple identifies the ANOVA component
        """
        # init, constant term + own term
        fi = {'const': ((-1) ** len(perm), self.__E),
              'var': [(1, perm)]}

        # add all lower order terms with alternating coefficient
        for k in xrange(len(perm) - 1):
            pperms = it.combinations(list(perm), r=k + 1)
            for pperm in pperms:
                fi['var'] += [((-1) ** (len(perm) - len(pperm)), pperm)]

        return fi

    def __decompose(self):
        """
        Computes the ANOVA components for the given function
        when the marginalized parts are available
        """
        fis = {}

        # add higher order terms
        for perm in self.__expectation_funcs.keys():
            fis[perm] = self.__combine_terms(perm)

        return fis

    def doDecomposition(self):
        """
        Computes the ANOVA decomposition for the given sparse grid function
        and the corresponding marginal distributions
        """
        if not self.__anova_components:
            # discretize variance function and use the same grid for all
            # following computations
            self.__discretize()
            # compute mean and variance
            self.__computeMean()
            self.__computeVariance()
            # marginalize
            self.__expectation_funcs = self.__expectation_functions()
            # combine marginalized terms
            self.__anova_components = self.__decompose()

    def __evalHigherOrderComponent(self, fi, x):
        """
        Evaluate the higher order components
        @param fi: linear combination of marginalized functions
        @param x: DataVector coordinate to be evaluated
        """
        # constant term
        sign, f0 = fi['const']
        s = sign * f0

        # higher order terms terms
        for sign, pperm in fi['var']:
            grid, alpha = self.__expectation_funcs[pperm]
            p = DataVector([x[ix] for ix in pperm])
            val = evalSGFunction(grid, alpha, p)
            s += sign * val

        return s

    def eval(self, x):
        """
        Evaluate the ANOVA decomposition at the given position
        @param x: coordinates to be evaluated
        """
        if not self.__anova_components:
            raise Exception('Do the decomposition first')

        # type check
        if isNumerical(x):
            x = np.array(x, dtype='float')

        # evaluation function

        # add constant term
        s = self.__E

        # add higher order terms
        for components in self.__anova_components.values():
            s += self.__evalHigherOrderComponent(components, x)

        return s

    def evalComponent(self, perm, x):
        """
        Evaluate a single ANOVA component
        @param perm: identifier
        @param x: coordinates
        """
        if len(perm) == 0:
            return self.__E
        elif perm in self.__anova_components:
            fi = self.__anova_components[perm]
            return self.__evalHigherOrderComponent(fi, x)
        else:
            raise AttributeError('The component %s does not exist' % (perm,))

    def getVarianceDecomposition(self):
        """
        Compute the variance of each ANOVA component
        """
        # check if this taks has alread been done
        if self.__variance_components:
            return self.__variance_components

        # do the anova decompisition first
        self.doDecomposition()

        # initialization
        vis = {}
        self.__variance_components = {}

        # run over all available permutations and compute the variance
        keys = self.__ap.keys()
        for perm in self.getSortedPermutations(self.__anova_components.keys()):
            # get the sparse grid function
            grid, alpha = self.__expectation_funcs[perm]

            # prepare parameter set
            dd = [keys[d] for d in perm]
            params = self.__ap.getSubset(dd)

            # transformation, pdf and volume
            T = params.getJointTransformation()
            U = params.getIndependentJointDistribution()

            # estimate variance
            def f(p, val):
                return val ** 2

            vi, err = self.__estimation.estimate(T.vol(), grid, alpha, f, U, T)
            vi -= self.__E ** 2
            # store the variance
            vis[perm] = vi

            if self._verbose:
                print "Estimated V[%s]: %g, L2 err = %g" % (perm, vi, err)
                print "-" * 60

            # add lower order components
            fi = self.__anova_components[perm]
            for sign, pperm in fi['var']:
                if len(pperm) < len(perm):
                    vi += sign * vis[pperm]

            self.__variance_components[perm] = vi

        return self.__variance_components

    def getSobolIndices(self):
        """
        Computes the relative influence of one single parameter combination
        to the whole variance.
        """
        # make sure that there exists a variance
        if self.__V > 0:
            vis = self.getVarianceDecomposition()
            ans = dict([(perm, vi / self.__V) for perm, vi in vis.items()])
        else:
            ans = {}
            for k in xrange(self.__nk):
                perms = it.combinations(range(self.__dim), r=k + 1)
                for perm in perms:
                    ans[perm] = 0.0
            if self.__has_highest_order_term:
                ans[tuple(range(self.__dim))] = 0.0

        return ans

    def getMainEffects(self):
        """
        Just compute the sobol indices without interactions
        """
        sobol = self.getSobolIndices()
        me = {}
        for perm in sobol.keys():
            if len(perm) == 1:
                me[perm] = sobol[perm]
        return me

    def getTotalEffects(self):
        """
        Compute the sum of all sobol indices where one single parameter
        is part of the interaction
        """
        sobol = self.getSobolIndices()
        te = {}
        for perm in sobol.keys():
            if len(perm) == 1:
                s = [sobol[pperm] for pperm in sobol.keys()
                     if perm[0] in pperm]
                te[perm] = sum(s)

        return te

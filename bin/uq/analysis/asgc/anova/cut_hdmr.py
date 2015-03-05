#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org
#
"""
@file    cut_hdmr.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 23 12:23:32 2013

@brief ANOVA-(H)igh (D)imensional (M)odel (R)epresentation
decomposition for a Sparse Grid

@version  0.1

"""

from pysgpp import (DataVector,
                    createOperationEval,
                    createOperationHierarchisation)

from bin.uq.quadrature import doMarginalize

import numpy as np
import itertools as it


class CutHDMR(object):
    """
    The CutHDMR class
    """

    def __init__(self, grid, alpha, U, trans, E, nk):
        self.__grid = grid
        self.__alpha = alpha
        self.__U = U
        self.__trans = trans
        self.__xlim = U.getBounds()
        self.__dim = len(U.getBounds())

        self.__E = E

        # check if highest order term is required
        self.__has_highest_order_term = nk == self.__dim

        if nk > self.__dim - 1:
            self.__nk = self.__dim - 1
        else:
            self.__nk = nk

        self.__expectation_funcs = None
        self.__anova_components = None
        self.__variance_components = None

        self.__verbose = True


    def __interpolate(self, grid, f):
        gs = grid.getStorage()

        # create new surpluss vector
        n_alpha = DataVector(gs.size())
        n_alpha.setAll(0.0)

        for i in xrange(gs.size()):
            gp = gs.get(i)
            p = [gp.abs(j) for j in xrange(gs.dim())]
            n_alpha[i] = f(p)

        createOperationHierarchisation(grid).doHierarchisation(n_alpha)

        return n_alpha


    # def __getSortedPermutations(self, keys):
    #     ans = [tuple()] * len(keys)
    #     ix = 0
    #     for x in sorted(np.unique(map(len, keys))):
    #         for ck in filter(lambda k: len(k) == x, keys):
    #             ans[ix] = ck
    #             ix += 1

    #     return ans


    def __expectation_functions(self):
        # create function to be interpolated
        opEval = createOperationEval(self.__grid)

        if self.__nk < 0 or self.__nk > self.__dim:
            raise AttributeError('The truncated order has to be in [0, %i]' % self.__dim)

        # init
        expec = {}

        # add higher order terms
        for k in range(self.__nk):
            perms = it.combinations(range(self.__dim), r=k + 1)
            for perm in perms:
                # select dimensions to be integrated
                dd = filter(lambda x: all([x != y for y in perm]),
                            xrange(self.__dim))

                if self.__verbose:
                    print "Explore %s, Integrate: %s" % (perm, dd)

                # -----------------------------------------------
                # Make sure that perm and dd are disjoint sets
                # covering the whole set
                assert sorted(list(perm) + dd) == range(self.__dim)
                assert len(dd) == self.__dim - len(perm)
                assert len(dd) > 0
                # -----------------------------------------------

                # @todo: do better interpolation considering adaptivity + distribution
                def f(p):
                    val = opEval.eval(self.__alpha, DataVector(p))
                    q = self.__U.pdf(self.__trans(p), marginal=True)
                    return val * np.prod(q[np.array(dd)])

                n_alpha = self.__interpolate(self.__grid, f)

                # do marginalization
                # => returns sparse grid function (grid, alpha)
                m_grid, m_alpha = doMarginalize(self.__grid, n_alpha, dd)

                # consider prefactor A = |Omega_i|
                A = float(np.prod(np.diff(self.__xlim[np.array(dd)])))
                m_alpha.mult(A)

                # store function
                expec[perm] = m_grid, m_alpha

        if self.__has_highest_order_term:
            expec[tuple(range(self.__dim))] = self.__grid, self.__alpha

        return expec


    def __combine_terms(self, perm):
        # init, constant term + own term
        fi = {'const': ((-1)**len(perm), self.__E),
              'var': [(1, perm)]}

        # add all lower order terms
        for k in xrange(len(perm) - 1):
            pperms = it.combinations(list(perm), r=k + 1)
            for pperm in pperms:
                fi['var'] += [((-1)**(len(perm) - len(pperm)), pperm)]

        return fi


    def __decompose(self):
        fis = {}

        # add higher order terms
        for perm in self.__expectation_funcs.keys():
            fis[perm] = self.__combine_terms(perm)

        return fis


    def doDecomposition(self):
        if not self.__anova_components:
            self.__expectation_funcs = self.__expectation_functions()
            self.__anova_components = self.__decompose()


    def __evalHigherOrderComponent(self, fi, x):
        # constant term
        sign, f0 = fi['const']
        s = sign * f0

        # higher order terms terms
        for sign, pperm in fi['var']:
            grid, alpha = self.__expectation_funcs[pperm]
            opEval = createOperationEval(grid)

            ixs = np.array(pperm, dtype='int')

            s += sign * opEval.eval(alpha, DataVector(x[ixs]))

        return s


    def eval(self, x):
        if not self.__anova_components:
            raise Exception('Do the decomposition first')

        # evaluation function

        # add constant term
        s = self.__E

        # add higher order terms
        for components in self.__anova_components.values():
            s += self.__evalHigherOrderComponent(components, x)

        return s


    def evalComponent(self, perm, x):
        if len(perm) == 0:
            return self.__E
        elif self.__anova_components.has_key(perm):
            fi = self.__anova_components[perm]
            return self.__evalHigherOrderComponent(fi, x)
        else:
            raise AttributeError('The component %s does not exist' % (perm, ))


    def getVarianceDecomposition(self, n=10000, *args, **kws):
        if self.__variance_components:
            return self.__variance_components
        else:
            self.doDecomposition(*args, **kws)

            xlim = np.array(self.__U.getBounds(), dtype='float')

            pos = np.random.rand(n, self.__dim)

            vis = {}

            for perm in self.__anova_components.keys():
                fi = self.__anova_components[perm]

                # do monte carlo integration

                # volume
                A = np.prod(np.diff(xlim[np.array(perm, dtype='int')]))

                # sample the domain
                s = 0.0
                for p in pos:
                    val = self.__evalHigherOrderComponent(fi, p)
                    q = self.__U.pdf(self.__trans(p), marginal=True)
                    s += val ** 2 * np.prod(q[np.array(perm)])

                # integral = ith variance component
                vi = A * s / len(pos)
                vis[perm] = vi


            self.__variance_components = vis

            return vis


    def getSobolIndices(self, V, *args, **kws):
        vis = self.getVarianceDecomposition(*args, **kws)
        ans = dict([(perm, vi / V) for perm, vi in vis.items()])
        return ans


    def getTotalEffects(self, *args, **kws):
        me = self.getSobolIndices(*args, **kws)
        te = {}
        for perm in me.keys():
            if len(perm) == 1:
                s = [me[pperm] for pperm in me.keys() if perm[0] in pperm]
                te[perm] = sum(s)

        return te


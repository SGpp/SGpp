# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# This file is part of sgpp, a program package making use of spatially
# adaptive sparse grids to solve numerical problems
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de),
#               2007-2009 Dirk Pflueger (pflueged@in.tum.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from pysgpp import *
from pysgpp.extensions.datadriven import *

def ApplyA(B, C, alpha, result, x, l):
    temp = DataVector(x.getNrows())
    M = x.getNrows();

    B.multBTrans(alpha, x, temp)
    B.multB(temp, x, result)

    temp = DataVector(len(alpha))

#    C.updown(alpha, temp)
#    result.axpy(M*l, temp)

    #for i in xrange(len(alpha)):
    #    temp[i] = alpha[i]

    result.axpy(M*l, alpha)


def sd(y, alpha, grid, x, imax, epsilon, l):
    B = OpB(grid)
    C = OpLaplace(grid)

    b = DataVector(len(alpha))
    B.multB(y, x, b)

    epsilon2 = epsilon*epsilon

    i = 1

    # r0 = b - Ax
    temp = DataVector(len(alpha))
    ApplyA(B, C, alpha, temp, x, l)
    r = DataVector(b)
    r.sub(temp)

    d = r.dotProduct(r)
    d0 = d * epsilon * epsilon
    print "delta_0 %g" % d0

    while i < imax+1 and d > d0:
        ApplyA(B, C, r, temp, x, l)
        a = d/r.dotProduct(temp)

        # x = x + ar
        alpha.axpy(a, r)

        if i % 100 == 0:
            ApplyA(B, C, alpha, temp, x, l)
            r.copyFrom(b)
            r.sub(temp)
        else:
            # r = r - aq
            r.axpy(-a, temp)

        d = r.dotProduct(r)
        print "delta: %g" % d
        i += 1
    print i
    print imax
    print d0
    print d


def cg(y, alpha, grid, x, imax, epsilon, l, verbose=True):
    B = OpB(grid)
    C = OpLaplace(grid)

    b = DataVector(len(alpha))
    B.multB(y, x, b)

    epsilon2 = epsilon*epsilon

    i = 1
    temp = DataVector(len(alpha))
    q = DataVector(len(alpha))

    ApplyA(B, C, alpha, temp, x, l)

    r = DataVector(b)
    r.sub(temp)

    d = DataVector(r)

    delta_old = 0.0
    delta_new = r.dotProduct(r)
    delta_0 = delta_new*epsilon2

    if verbose:
        print "delta_0 %g" % delta_0

    while (i < imax+1) and (delta_new > delta_0):
    	# q = A*d
        ApplyA(B, C, d, q, x, l)
	    # a = d_new / d.q
        a = delta_new/d.dotProduct(q)

        # x = x + a*d
        alpha.axpy(a, d)

        if i % 50 == 0:
	    # r = b - A*x
            ApplyA(B, C, alpha, temp, x, l)
            r.copyFrom(b)
            r.sub(temp)
        else:
            # r = r - a*q
            r.axpy(-a, q)

        if verbose:
            print "delta: %g" % delta_new

        delta_old = delta_new
        delta_new = r.dotProduct(r)
        beta = delta_new/delta_old

        d.mult(beta)
        d.add(r)

        i += 1

    if verbose:
        print i
        print imax
        print delta_0
        print delta_new


def BiCGStab(b, alpha, imax, epsilon, ApplyMatrix, verbose=True):
    # nach http://www.iue.tuwien.ac.at/phd/heinreichsberger/node70.html
    # http://www.numerik.math.tu-graz.ac.at/kurse/lgs/SIMNET6.pdf

    i = 1
    epsilon2 = epsilon*epsilon
    #temp = DataVector(len(alpha))

    # Choose x0
    alpha.setAll(0.0)

    # Calculate r0
    r = DataVector(len(alpha))
    ApplyMatrix(alpha, r)
    r.sub(b)

    delta_0 = r.dotProduct(r)*epsilon2
    delta = 0.0
    if verbose:
        print "delta_0 %g" % delta_0

    # Choose r0 as r
    r0 = DataVector(r)
    # Set p as r0
    p = DataVector(r0)

    rho = r0.dotProduct(r)

    s = DataVector(len(alpha))

    while i < imax:
        # s  = Ap
        ApplyMatrix(p, s)

        sigma = s.dotProduct(r0)
        if sigma == 0.0:
            break

        a = rho/sigma

        #w = r - a*s
        w = DataVector(r)
        w.axpy(-a, s)

        #v = Aw
        v = DataVector(len(alpha))
        ApplyMatrix(w, v)

        omega = v.dotProduct(w) / v.dotProduct(v)

        #x = x - a*p - omega*w
        alpha.axpy(-a, p)
        alpha.axpy(-omega, w)

        #r = r - a*s - omega*v
        r.axpy(-a, s)
        r.axpy(-omega, v)

        rho_new = r.dotProduct(r0)

        delta = r.dotProduct(r)
        if verbose:
            print "delta: %g" % delta
        # stoppe falls genauigkeit
        if delta < delta_0:
            break

        beta = rho_new/rho * a/omega
        rho = rho_new

        # p = r + beta*(p - omega*s)
        p.axpy(-omega, s)
        p.mult(beta)
        p.add(r)

    return (i, delta)


#-------------------------------------------------------------------------------
## Conjugated Gradient method for sparse grids, solving A.alpha=b.
# The resulting vector is stored in alpha.
# @param b RHS of equation
# @param alpha vector of unknowns
# @param imax max. number of iterations (abort, if reached)
# @param epsilon accuracy requirements (reduce initial norm of residuum |delta_0|
#        below epsilon*|delta_0|)
# @param ApplyMatrix procedure that applies A to a vector
# @param reuse starting vector is 0 by default. If true, use current values in alpha
# @param verbose verbose output (default False)
# @param max_threshold maximal threshold
# @return tuple (number of iterations, final norm of residuum)
def cg_new(b, alpha, imax, epsilon, ApplyMatrix, reuse = False, verbose=True, max_threshold=None):
    if verbose:
        print "Starting Conjugated Gradients"

#    Apply B to y
#    b = DataVector(len(alpha))
#    B.multB(y, x, b)

    epsilon2 = epsilon*epsilon

    i = 0
    temp = DataVector(len(alpha))
    q = DataVector(len(alpha))
    delta_0 = 0.0

    # calculate residuum
    if reuse:
        q.setAll(0)
        ApplyMatrix(q, temp)
        r = DataVector(b)
        r.sub(temp)
        delta_0 = r.dotProduct(r)*epsilon2
    else:
        alpha.setAll(0)

    ApplyMatrix(alpha, temp)
    r = DataVector(b)
    r.sub(temp)

    # delta
    d = DataVector(r)

    delta_old = 0.0
    delta_new = r.dotProduct(r)

    if not reuse:
        delta_0 = delta_new*epsilon2

    if verbose:
        print "Starting norm of residuum: %g" % (delta_0/epsilon2)
        print "Target norm:               %g" % (delta_0)

    while (i < imax) and (delta_new > delta_0) and (max_threshold == None or delta_new > max_threshold):
        # q = A*d
        ApplyMatrix(d, q)
        # a = d_new / d.q
        a = delta_new/d.dotProduct(q)

        # x = x + a*d
        alpha.axpy(a, d)

        if i % 50 == 0:
        # r = b - A*x
            ApplyMatrix(alpha, temp)
            r.copyFrom(b)
            r.sub(temp)
        else:
            # r = r - a*q
            r.axpy(-a, q)

        delta_old = delta_new
        delta_new = r.dotProduct(r)
        beta = delta_new/delta_old

        if verbose:
            print "delta: %g" % delta_new

        d.mult(beta)
        d.add(r)

        i += 1

    if verbose:
        print "Number of iterations: %d (max. %d)" % (i, imax)
        print "Final norm of residuum: %g" % delta_new

    return (i,delta_new)

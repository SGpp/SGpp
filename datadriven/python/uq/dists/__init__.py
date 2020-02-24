# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
Statistical Distributions
===============================

This module contains univariate distribution functions

Supported univariate distributions
------------------------------------
normal
uniform
beta
lognormal
truncated normal
sparse grid denstiy estimation distribution
data dist

Supported multivariate distributions
------------------------------------
J (independent multivariate-distribution)
Corr (correlated multivariate-distribution)

"""

__version__ = "1.0"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from pysgpp.extensions.datadriven.uq.dists.Beta import Beta
from pysgpp.extensions.datadriven.uq.dists.Corr import Corr
from pysgpp.extensions.datadriven.uq.dists.CorrBeta import CorrBeta
from pysgpp.extensions.datadriven.uq.dists.DataDist import DataDist
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
from pysgpp.extensions.datadriven.uq.dists.J import J
from pysgpp.extensions.datadriven.uq.dists.Lognormal import Lognormal
from pysgpp.extensions.datadriven.uq.dists.TLognormal import TLognormal
from pysgpp.extensions.datadriven.uq.dists.TNormal import TNormal
from pysgpp.extensions.datadriven.uq.dists.Normal import Normal
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform
from pysgpp.extensions.datadriven.uq.dists.MultivariateNormal import MultivariateNormal

from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
from pysgpp.extensions.datadriven.uq.dists.LibAGFDist import LibAGFDist
from pysgpp.extensions.datadriven.uq.dists.KDEDist import KDEDist
try:
    from pysgpp.extensions.datadriven.uq.dists.NatafDist import NatafDist
    from pysgpp.extensions.datadriven.uq.dists.DTreesDist import DTreesDist
except:
    pass

# import optimization

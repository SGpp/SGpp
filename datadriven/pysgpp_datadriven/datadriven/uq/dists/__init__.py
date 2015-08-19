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

from Beta import Beta
from Corr import Corr
from CorrBeta import CorrBeta
from DataDist import DataDist
from Dist import Dist
from J import J
from Lognormal import Lognormal
from TNormal import TNormal
from Normal import Normal
from Uniform import Uniform
from MultivariateNormal import MultivariateNormal

from SGDEdist import SGDEdist
from LibAGFDist import LibAGFDist
from DTreesDist import DTreesDist
from GaussianKDEDist import GaussianKDEDist

import optimization

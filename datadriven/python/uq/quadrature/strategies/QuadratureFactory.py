# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.dists import TNormal, Normal, Uniform, SGDEdist
from pysgpp.extensions.datadriven.uq.quadrature.strategies.GaussHermiteQuadrature import GaussHermiteQuadrature
from pysgpp.extensions.datadriven.uq.quadrature.strategies.GaussLegendreQuadrature import GaussLegendreQuadrature
import warnings


class QuadratureFactory(object):

    @classmethod
    def findQuadratureStrategyByMeasure(cls, U):
        if U is None:
            return GaussLegendreQuadrature()
        if all([isinstance(Ui, Uniform) or isinstance(Ui, SGDEdist) 
                for Ui in U]):
            return GaussLegendreQuadrature()
#         elif all([isinstance(Ui, TNormal) or
#                   isinstance(Ui, Normal) for Ui in U]):
#             return GaussHermiteQuadrature()
        else:
            warnings.warn("no suitable quadrature strategy found. using gauss legendre")
            return GaussLegendreQuadrature()

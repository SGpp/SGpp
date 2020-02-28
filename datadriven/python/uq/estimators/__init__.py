# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
Estimators
==========================================

This module contains the different estimators for statistical moment
extraction of sparse grid functions

"""

__version__ = "0.1"

__all__ = []

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

# from CollocationPointsStrategy import CollocationPointsStrategy
# from SparseGridEstimationStrategy import SparseGridEstimationStrategy
# from SparseGridEstimator import SparseGridEstimator
# from MonteCarloStrategy import MonteCarloStrategy
# from PiecewiseConstantIntegralStrategy import PiecewiseConstantIntegralStrategy
# from DiscreteIntegralStrategy import DiscreteIntegralStrategy
# from AnalyticEstimationStrategy import AnalyticEstimationStrategy
# from MarginalAnalyticEstimationStrategy import MarginalAnalyticEstimationStrategy


from pysgpp.extensions.datadriven.uq.estimators.MCEstimator import MCEstimator
from pysgpp.extensions.datadriven.uq.estimators.AnalyticEstimationStrategy import AnalyticEstimationStrategy
from pysgpp.extensions.datadriven.uq.estimators.CollocationPointsStrategy import CollocationPointsStrategy
from pysgpp.extensions.datadriven.uq.estimators.MarginalAnalyticEstimationStrategy import MarginalAnalyticEstimationStrategy
from pysgpp.extensions.datadriven.uq.estimators.MonteCarloStrategy import MonteCarloStrategy
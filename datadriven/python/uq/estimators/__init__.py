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


from MCEstimator import MCEstimator
from AnalyticEstimationStrategy import AnalyticEstimationStrategy
from CollocationPointsStrategy import CollocationPointsStrategy
from MarginalAnalyticEstimationStrategy import MarginalAnalyticEstimationStrategy
from MonteCarloStrategy import MonteCarloStrategy
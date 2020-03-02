# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
Learner
==========================================

"""

__version__ = "0.1"

__all__ = ["Interpolant", "ANOVAInterpolant"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"

from pysgpp.extensions.datadriven.uq.learner.builder.CGSolverDescriptor import CGSolverDescriptor
from pysgpp.extensions.datadriven.uq.learner.builder.GridDescriptor import GridDescriptor
from pysgpp.extensions.datadriven.uq.learner.builder.RegressorSpecificationDescriptor import RegressorSpecificationDescriptor
# from pysgpp.extensions.datadriven.uq.learner.builder.SimulationLearnerBuilder import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.uq.learner.builder.StopPolicyDescriptor import StopPolicyDescriptor

from pysgpp.extensions.datadriven.uq.learner.Interpolant import Interpolant
from pysgpp.extensions.datadriven.uq.learner.Learner import Learner, LearnerEvents
from pysgpp.extensions.datadriven.uq.learner.Regressor import Regressor
# from SimulationLearner import SimulationLearner

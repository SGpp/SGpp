# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
Parameter Specification
==========================================

"""

__version__ = "0.1"

__all__ = ["ParameterBuilder", "Parameter", "ParameterSet",
           "UncertainParameterBuilder", "UncertainParameter"
           "DeterministicParameterBuilder", "DeterministicParameter"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"


from pysgpp.extensions.datadriven.uq.parameters.DeterministicParameter import DeterministicParameter
from pysgpp.extensions.datadriven.uq.parameters.Parameter import Parameter
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import (ParameterBuilder,
                              UncertainParameterBuilder,
                              DeterministicParameterBuilder)
from pysgpp.extensions.datadriven.uq.parameters.ParameterSet import ParameterSet
from pysgpp.extensions.datadriven.uq.parameters.UncertainParameter import UncertainParameter

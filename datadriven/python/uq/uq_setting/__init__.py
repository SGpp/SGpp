# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

"""
UQ Setting
===============================

This module contains a framework for non-intrusive UQ. The interface
is realized using a fluid-interface-pattern. A UQSetting can be built
using the

UQBuilder

which has the method andGetResult() returning a

UQSetting

object, which is capable of restoring preprocessing, simulation and
postprocessing results.
"""

__version__ = "0.1"

__all__ = ["UQBuilder", "UQSetting",
           "UQSettingFormatter", "UQSpecification"]

__author__ = "Fabian Franzelin, fabian.franzelin@ipvs.uni-stuttgart.de"


from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder
from pysgpp.extensions.datadriven.uq.uq_setting.UQSetting import UQSetting, UQSampleType
from pysgpp.extensions.datadriven.uq.uq_setting.UQSettingAdapter import UQSettingAdapter
from pysgpp.extensions.datadriven.uq.uq_setting.UQSettingFormatter import UQSettingFormatter
from pysgpp.extensions.datadriven.uq.uq_setting.UQSpecification import UQSpecification
from pysgpp.extensions.datadriven.uq.uq_setting.UQSettingTools import findEquivalent


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


from UQBuilder import UQBuilder
from UQSetting import UQSetting, UQSampleType
from UQSettingAdapter import UQSettingAdapter
from UQSettingFormatter import UQSettingFormatter
from UQSpecification import UQSpecification
from UQSettingTools import findEquivalent


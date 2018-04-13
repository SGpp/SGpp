# !/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    ASGC.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 23 12:58:31 2013

@brief   Adaptive Sparse Grid Collocation method for UQ

@version  0.1

"""
from pysgpp.extensions.datadriven.uq.estimators import MCEstimator
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp import (DataVector, DataMatrix)
import numpy as np

from pysgpp.extensions.datadriven.uq.analysis import Analysis


class MCAnalysis(Analysis):
    """
    The MCAnalysis class
    """

    def __init__(self, params, samples, estimator=None,
                 npaths=100):
        """
        Constructor
        @param params: ParameterSet
        @param samples: dictionary {<time step>: {<Sample>: value}}
        """
        Analysis.__init__(self, ts=samples.keys())
        self.__params = params
        self.__samples = samples
        if estimator is None:
            self.__estimator = MCEstimator(npaths)
        else:
            self.__estimator = estimator

    def computeMean(self, iteration, qoi, t):
        # do the computation
        values = np.array(self.__samples[t].values())
        return self.__estimator.mean(values)

    def computeVar(self, iteration, qoi, t):
        values = np.array(self.__samples[t].values())
        return self.__estimator.var(values)

# -----------------------------------------------------------------------------

    def computeMoments(self, ts=None):
        names = ['time',
                 'iteration',
                 'grid_size',
                 'mean',
                 'mean_err',
                 'meanConfidenceIntervalBootstrapping_lower',
                 'meanConfidenceIntervalBootstrapping_upper',
                 'var',
                 'var_err',
                 'varConfidenceIntervalBootstrapping_lower',
                 'varConfidenceIntervalBootstrapping_upper']
        # parameters
        ts = self.__samples.keys()
        nrows = len(ts)
        ncols = len(names)
        data = DataMatrix(nrows, ncols)
        v = DataVector(ncols)

        row = 0
        for t in np.sort(ts):
            v.setAll(0.0)
            mean = self.mean(ts=[t], iterations=[0])
            var = self.var(ts=[t], iterations=[0])
            numSamples = len(self.__samples[t].values())

            v[0] = t
            v[1] = 0
            v[2] = numSamples
            v[3], v[4] = mean["value"], mean["err"]
            v[5], v[6] = mean["confidence_interval"]
            v[7], v[8] = var["value"], var["err"]
            v[8], v[9] = var["confidence_interval"]

            # write results to matrix
            data.setRow(row, v)
            row += 1

        return {'data': data,
                'names': names}

    def writeMoments(self, filename):
        stats = self.computeMoments()
        stats['filename'] = filename + ".moments.arff"
        writeDataARFF(stats)

# -----------------------------------------------------------------------------

    def estimateDensity(self, ts=[0], dtype="kde", config={}):
        if len(ts) == 1:
            return self._estimateDensityByConfig(dtype, self.__samples[ts[0]].values(), config)

        ans = {}
        for t, values in time_dependent_values.items():
            ans[t] = self._estimateDensityByConfig(dtype, values, config)

        return ans

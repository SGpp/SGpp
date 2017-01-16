# !/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
from pysgpp.extensions.datadriven.uq.dists.KDEDist import KDEDist
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
from scipy.stats.mstats_basic import moment
"""
@file    ASGC.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 23 12:58:31 2013

@brief   Adaptive Sparse Grid Collocation method for UQ

@version  0.1

"""
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
import numpy as np

class AnalysisHash(object):

    def __init__(self):
        self._moments = {}

    def reset(self):
        self._moments = {}

    def getKey(self, iteration, qoi, t, idd):
        return (iteration, qoi, t, idd)

    def hasMoment(self, iteration, qoi, t, idd='mean'):
        return self.getKey(iteration, qoi, t, idd) in self._moments

    def setMoment(self, iteration, qoi, t, idd, value):
        key = self.getKey(iteration, qoi, t, idd)
        self._moments[key] = value

    def getMoment(self, iteration, qoi, t, idd):
        key = self.getKey(iteration, qoi, t, idd)
        if self.hasMoment(iteration, qoi, t, idd):
            return self._moments[key]
        else:
            return None


class Analysis(object):
    """
    the analysis class
    """

    def __init__(self, qoi="_", ts=[0], iterations=[0]):
        # initialize dicts
        self._qoi = qoi
        self._ts = ts
        self._iterations = iterations
        self._moments = AnalysisHash()
        self._verbose = True

    def setVerbose(self, verbose):
        self._verbose = verbose

# -----------------------------------------------------------------------------

    def computeMean(self, iteration, qoi, t):
        raise NotImplementedError()

    def mean(self, iterations=None, ts=None):
        """
        compute means
        @return: dictionary, {<iteration>: {<time>: (mean, err)}}
        """
        if iterations is None:
            iterations = self._iterations

        if ts is None:
            ts = self._ts

        ans = {}
        for iteration in iterations:
            if len(ts) > 1:
                ans[iteration] = {}
            for i, t in enumerate(ts):
                # compute mean
                if self._verbose:
                    print "-" * 60
                    print "Estimate E[t = %g] (%i/%i), iteration = %s:" % \
                        (t, i + 1, len(self._ts), iteration),

                if not self._moments.hasMoment(iteration, self._qoi, t, 'mean'):
                    moment = self.computeMean(iteration, self._qoi, t)
                    self._moments.setMoment(iteration, self._qoi, t,
                                            'mean', moment)
                else:
                    moment = self._moments.getMoment(iteration, self._qoi,
                                                     t, 'mean')

                if self._verbose:
                    print "value = %g (err=%g)" % moment

                if len(ts) > 1:
                    ans[iteration][t] = moment
                else:
                    ans[iteration] = moment

        # remove dict structure if there are just one element
        if len(iterations) == 1:
            ans = ans[iterations[0]]

        return ans

# -----------------------------------------------------------------------------

    def computeVar(self, iteration, qoi, t):
        raise NotImplementedError()

    def var(self, iterations=None, ts=None):
        """
        Compute the variance
        @return: dictionary, {<iteration>: {<time>: variance}}
        """
        if iterations is None:
            iterations = self._iterations

        if ts is None:
            ts = self._ts

        ans = {}
        for iteration in iterations:
            if len(ts) > 1:
                ans[iteration] = {}
            for i, t in enumerate(ts):
                # compute variance
                if self._verbose:
                    print "-" * 60
                    print "Estimate V[t = %g] (%i/%i), iteration = %s:" % \
                        (t, i + 1, len(self._ts), iteration),

                if not self._moments.hasMoment(iteration, self._qoi, t, 'var'):
                    moment = self.computeVar(iteration, self._qoi, t)
                    self._moments.setMoment(iteration, self._qoi, t,
                                            'var', moment)
                else:
                    moment = self._moments.getMoment(iteration, self._qoi,
                                                     t, 'var')

                if self._verbose:
                    print "value = %g (err=%g)" % moment

                if len(ts) > 1:
                    ans[iteration][t] = moment
                else:
                    ans[iteration] = moment

        # remove dict structure if there are just one element
        if len(iterations) == 1:
            ans = ans[iterations[0]]

        return ans

# -----------------------------------------------------------------------------

    def confidenceInterval(self, iterations=None, ts=None):
        """
        Compute the variance
        @return: dictionary, {<iteration>: {<time>: variance}}
        """
        if iterations is None:
            iterations = self._iterations

        if ts is None:
            ts = self._ts

        ans = {}
        for iteration in iterations:
            if len(ts) > 1:
                ans[iteration] = {}
            for i, t in enumerate(self._ts):
                # compute variance
                if self._verbose:
                    print "-" * 60
                    print "Estimate confidence interval for t = %g (%i/%i), iteration = %s" % \
                        (t, i + 1, len(self._ts), iteration)

                confidenceInterval = self.computeConfidenceInterval(iteration, self._qoi, t)

                if len(ts) > 1:
                    ans[iteration][t] = confidenceInterval
                else:
                    ans[iteration] = confidenceInterval

        # remove dict structure if there are just one element
        if len(iterations) == 1:
            ans = ans[iterations[0]]

        return ans
        
        
    def computeConfidenceInterval(self, iteration, qoi, t):
        raise NotImplementedError()

# -----------------------------------------------------------------------------

    def computeMoments(self, iterations=None, ts=None):
        raise NotImplementedError()

    def writeMoments(self, filename, *args, **kws):
        stats = self.computeMoments(*args, **kws)
        stats['filename'] = filename + ".moments.arff"
        writeDataARFF(stats)

# -----------------------------------------------------------------------------

    def _estimateDensityByConfig(self, dtype, samples, config={}):
        if dtype == "kde":
            # compute bounds of samples
            bounds = np.ndarray((1, 2))
            bounds[:, 0] = np.min(samples)
            bounds[:, 1] = np.max(samples)
            return KDEDist(samples, bounds=bounds)
        elif dtype == "sgde":
            # compute bounds of samples
            bounds = np.ndarray((1, 2))
            bounds[:, 0] = np.min(samples)
            bounds[:, 1] = np.max(samples)
            return SGDEdist.byLearnerSGDEConfig(samples, bounds, config=config)
        else:
            raise AttributeError("density estimation type %s is not known. Select one in [gaussianKDE, sgde]")

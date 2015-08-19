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
from pysgpp_datadriven.uq.tools import writeDataARFF


class AnalysisHash(object):

    def __init__(self):
        self._moments = {}

    def reset(self):
        self._moments = {}

    def hasMoment(self, iteration, qoi, t, idd='mean'):
        return iteration in self._moments and \
            qoi in self._moments[iteration] and \
            t in self._moments[iteration][qoi] and \
            idd in self._moments[iteration][qoi][t]

    def setMoment(self, iteration, qoi, t, idd, value):
        # check if dictionary is available
        if iteration not in self._moments:
            self._moments[iteration] = {}
        if qoi not in self._moments[iteration]:
            self._moments[iteration][qoi] = {}
        if t not in self._moments[iteration][qoi]:
            self._moments[iteration][qoi][t] = {}

        # store the value in the dictionary
        self._moments[iteration][qoi][t][idd] = value

    def getMoment(self, iteration, qoi, t, idd):
        if self.hasMoment(iteration, qoi, t, idd):
            return self._moments[iteration][qoi][t][idd]
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
        @return: dictionary, means
        """
        if iterations is None:
            iterations = self._iterations
        if ts is None:
            ts = self._ts

        ans = {}
        for iteration in iterations:
            ans[iteration] = {}
            for i, t in enumerate(self._ts):
                # compute mean
                if self._verbose:
                    print "-" * 60
                    print "Estimate E[t = %g] (%i/%i), iteration = %s" % \
                        (t, i + 1, len(self._ts), iteration)

                if not self._moments.hasMoment(iteration, self._qoi, t, 'mean'):
                    moment = self.computeMean(iteration, self._qoi, t)
                    self._moments.setMoment(iteration, self._qoi, t,
                                            'mean', moment)
                else:
                    moment = self._moments.getMoment(iteration, self._qoi,
                                                     t, 'mean')

                ans[iteration][t] = moment

        # remove dict structure if there are just one element
        if len(iterations) == 1:
            ans = ans[iterations[0]]
        if len(ts) == 1:
            ans = ans[ts[0]]

        return ans

# -----------------------------------------------------------------------------

    def computeVar(self, iteration, qoi, t):
        raise NotImplementedError()

    def var(self, iterations=None, ts=None):
        """
        Compute the variance
        @return: dictionary, variances
        """
        if iterations is None:
            iterations = self._iterations
        if ts is None:
            ts = self._ts

        ans = {}
        for iteration in iterations:
            ans[iteration] = {}
            for i, t in enumerate(self._ts):
                # compute variance
                if self._verbose:
                    print "-" * 60
                    print "Estimate V[t = %g] (%i/%i), iteration = %s" % \
                        (t, i + 1, len(self._ts), iteration)

                if not self._moments.hasMoment(iteration, self._qoi, t, 'var'):
                    moment = self.computeVar(iteration, self._qoi, t)
                    self._moments.setMoment(iteration, self._qoi, t,
                                            'var', moment)
                else:
                    moment = self._moments.getMoment(iteration, self._qoi,
                                                     t, 'var')

                ans[iteration][t] = moment

        # remove dict structure if there are just one element
        if len(iterations) == 1:
            ans = ans[iterations[0]]
        if len(ts) == 1:
            ans = ans[ts[0]]

        return ans

# -----------------------------------------------------------------------------

    def computeMoments(self, iterations=None, ts=None):
        raise NotImplementedError()

    def writeMoments(self, filename, *args, **kws):
        stats = self.computeMoments(*args, **kws)
        stats['filename'] = filename + ".moments.arff"
        writeDataARFF(stats)

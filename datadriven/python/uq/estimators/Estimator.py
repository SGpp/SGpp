# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


class Estimator(object):

    def mean(self, *args, **kws):
        """
        Estimate the mean
        @return: tuple(moment, error)
        """
        raise NotImplementedError()

    def var(self, *args, **kws):
        """
        Estimate the variance
        """
        raise NotImplementedError()

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder

class Model(object):

    def __init__(self, params, filename=None):
        self.num_results = 1
        self.params = params
        self.filename = filename

        # load uq setting
        self.uqBuilder = UQBuilder()
        if filename is not None:
            self.uqBuilder.fromFile(filename)

        self.buildUQSetting()
        self.uqSetting = self.uqBuilder.andGetResult()

    def buildUQSetting(self):
        raise NotImplementedError

    def evaluate(self, samples, distributed=False):
        self.uqSetting.runSamples(samples, dist=distributed)

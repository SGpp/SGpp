
'''
Created on Feb 6, 2015

@author: franzefn
'''

from interpolationAlgorithm import InterpolationAlgorithm
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
import os
import tempfile, uuid, json


class EstimateDensityAlgorithm(InterpolationAlgorithm):

    def __init__(self, trainSamples, config={}):
        self.config = config.copy()
        self.trainSamples = trainSamples

        if "grid_filename" in self.config:
            self.gridfile = self.config["grid_filename"]
        else:
            # get temp directory
            self.gridfile = os.path.join(tempfile.gettempdir(),
                                         "sgde-grid-%s.grid" % str(uuid.uuid4()))
            self.config["grid_filename"] = self.gridfile

        # remove refinement if there is some
        if "refinement_numSteps" in self.config:
            self.config["refinement_numSteps"] = 0

    def computeHierarchicalCoefficients(self, grid, alpha, addedGridPoints=None):
        # serialize grid
        fd = open(self.gridfile, "w")
        fd.write(grid.serialize())
        fd.close()

        # estimate the density on the given grid
        dist = SGDEdist.byLearnerSGDEConfig(self.trainSamples, config=self.config)

        return dist.alpha

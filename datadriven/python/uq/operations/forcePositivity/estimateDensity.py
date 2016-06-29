
'''
Created on Feb 6, 2015

@author: franzefn
'''

import os
import numpy as np

from interpolationAlgorithm import InterpolationAlgorithm
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
import tempfile, uuid, json
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import dehierarchize, hierarchize
from pysgpp import (DataVector,
                    DataMatrix,
                    createOperationLTwoDotExplicit,
                    RegularizationType_Identity,
                    RegularizationType_Laplace,
                    createOperationIdentity,
                    createOperationLaplace,
                    createOperationMultipleEval)
from pysgpp.pysgpp_swig import createOperationLTwoDotProduct
from pysgpp.extensions.datadriven.tools import writeGrid


class EstimateDensityAlgorithm(InterpolationAlgorithm):

    def __init__(self, trainSamples):
        self.trainSamples = DataMatrix(trainSamples)
        self.lmbd = lmbd
        self.regularizationType = regularizationType

    def computeHierarchicalCoefficients(self, grid, alpha, addedGridPoints=None):
        nodalValues = dehierarchize(grid, alpha)
        ixs = np.array([], dtype="int")
        for i, yi in enumerate(nodalValues):
            if yi < 0:
                ixs = np.append(ixs, i)
                nodalValues[i] = 0.0

        if len(ixs) > 0:
            # compute the coefficients for each grid point by estimating
            # the local density
            config = {'grid_filename': "positive.grid",
                      "regularization_type": "Laplace",
                      "crossValidation_enable": True,
                      "crossValidation_enable": True,
                      "crossValidation_kfold": 5,
                      "crossValidation_silent": True}

            writeGrid(config['grid_filename'], grid)
            alpha = SGDEdist.byLearnerSGDEConfig(self.trainSamples.array(), bounds=None, config=config).alpha.array()

            # check if the coefficients of the new grid points are positive
            if addedGridPoints is not None:
                gs = grid.getStorage()
                assert all([alpha[gs.getSequenceNumber(gp)] > -1e-13 for gp in addedGridPoints])
        return alpha

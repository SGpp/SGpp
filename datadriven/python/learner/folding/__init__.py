# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

##############################################################################
                                    #
#############################################################################

from pysgpp.extensions.datadriven.learner.folding.FoldingPolicy import FoldingPolicy
from pysgpp.extensions.datadriven.learner.folding.RandomFoldingPolicy import RandomFoldingPolicy
from pysgpp.extensions.datadriven.learner.folding.SequentialFoldingPolicy import SequentialFoldingPolicy
from pysgpp.extensions.datadriven.learner.folding.StratifiedFoldingPolicy import StratifiedFoldingPolicy
from pysgpp.extensions.datadriven.learner.folding.FilesFoldingPolicy import FilesFoldingPolicy

__all__ = ['FilesFoldingPolicy',
            'FoldingPolicy',
            'RandomFoldingPolicy',
            'SequentialFoldingPolicy',
            'StratifiedFoldingPolicy',]
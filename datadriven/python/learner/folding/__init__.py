# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##############################################################################
                                    #
#############################################################################

from FoldingPolicy import FoldingPolicy
from RandomFoldingPolicy import RandomFoldingPolicy
from SequentialFoldingPolicy import SequentialFoldingPolicy
from StratifiedFoldingPolicy import StratifiedFoldingPolicy
from FilesFoldingPolicy import FilesFoldingPolicy

__all__ = ['FilesFoldingPolicy',
            'FoldingPolicy',
            'RandomFoldingPolicy',
            'SequentialFoldingPolicy',
            'StratifiedFoldingPolicy',]
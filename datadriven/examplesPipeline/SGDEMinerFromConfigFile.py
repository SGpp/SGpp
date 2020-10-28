#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## \page example_SGDEMinerFromConfigFile_py SGDE Miner
## This example demonstrates how a pysgpp.datadriven.SparseGridMiner
## is constructed using a configuration file and how it is then used

import sys
from pysgpp import DensityEstimationMinerFactory
from pysgpp import SparseGridMiner

# check for config file command line argument
if len(sys.argv) != 2:
	print("No or bad path given, aborting")
	sys.exit(1)

path = sys.argv[1]

# Instantiate the factory class
factory = DensityEstimationMinerFactory()

# Create miner using the factory and the given config file
miner = factory.buildMiner(path)

# Once configured, the miner object can now perform the learning
miner.learn(True)

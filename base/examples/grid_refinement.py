# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python
# @file grid_refinement.py
# This is an example of how to resume the learning process from a checkpoint
# without recalculating the last iteration 

# include all necessary classes
from bin.learner.LearnerBuilder import LearnerBuilder
from bin.controller import CheckpointController
from bin.controller import InfoToScreenRegressor, InfoToFile

if __name__=='__main__':
	# code automatically generated using TerminalController
	code  = """builder = LearnerBuilder()
builder.buildRegressor()
builder.withTrainingDataFromARFFFile('/path/to/dataset_train.arff.gz')
builder.withTestingDataFromARFFFile('/path/to/dataset_test.arff.gz')
builder = builder.withGrid()
builder.fromFile('/path/to/grid_file.grid.gz')
builder.withStartingIterationNumber(9)
builder = builder.withSpecification()
builder.withAdaptPoints(400)
builder.withLambda(0.00010000000000000000)
builder.withIdentityOperator()
builder = builder.withStopPolicy()
builder.withGridSizeLimit(50000)
builder = builder.withCGSolver()  
builder.withAccuracy(0.000100)
builder.withImax(500)
builder.withProgressPresenter(InfoToScreenRegressor())
builder.withProgressPresenter(InfoToFile('/path/to/stats_file.stats.gz'))
checkpointController = CheckpointController('checkpoint_file.checkpoint.gz', 
'/path/to')
builder.withCheckpointController(checkpointController)"""
	# create the Learner object
	exec code
	learner = builder.andGetResult()
	# update the values of errors per grid point
	dataset = learner.dataContainer
	trainSubset = dataset.getTrainDataset()
	testSubset = dataset.getTestDataset()
	learner.updateResults(learner.alpha, trainSubset, testSubset)
	# refine the grid
	learner.refineGrid()
	# resume the learning process
	learner.learnDataWithTest()
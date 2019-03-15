#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## @page example_resumingLearningProcess_py resumingLearningProcess.py
##
## This is an example of how to resume the learning process from a checkpoint
## without recalculating the last iteration.
##
## Sometimes one needs to resume the learning process using the saved grid and alpha
## files in order to refine the learner or to use some different learning
## parameters. And while the definition of the setup with command line arguments may
## be comfortable for the most common tasks, the default behaviour in this case
## will result in the  recalculating of the last learning iteration. In order to
## change this behaviour the direct interaction with Learner object is needed.
##
## This example shows how to create the Learner object using your command line
## arguments as well as how to omit the recomputation of the last iteration resuming
## the learning.
##
## In order to create a Learner object one should use the
## @link python.learner.LearnerBuilder.LearnerBuilder LearnerBuilder @endlink.
## You can let the
## @link python.controller.TerminalController.TerminalController TerminalController @endlink
## to generate the code from the command line arguments.
##
## Let's suppose you learning setup was defined by the following command line
## arguments
## @verbatim
## python $SGPP_PATH/misc/python/classifier.py -v --mode test --data \
## /path/to/dataset_train.arff.gz --test /path/to/dataset_test.arff.gz --border \
## --level 2 -r 0.0001 -C identity --adapt_points 400 --grid_limit 50000 \
## --regression --checkpoint /path/to/checkpoint_file.checkpoint.gz --stats \
## /path/to/stats_file.stats.gz --grid /path/to/grid_file.grid.gz --adapt_start 9 \
## -L 0.0001
## @endverbatim
##
## With @link python.controller.TerminalController.TerminalController TerminalController @endlink
## class you now can automatically generate the code to
## create the Learner object using LearnerBuilder. Just call the
## %TerminalController.py with all the command line arguments as above and an
## additional argument @c --generate
##
## @verbatim
## python $SGPP_PATH/datadriven/python/controller/TerminalController.py --generate -v --mode \
## test --data /path/to/dataset_train.arff.gz --test /path/to/dataset_test.arff.gz \
## --border --level 2 -r 0.0001 -C identity --adapt_points 400 --grid_limit 50000 \
## --regression --checkpoint /path/to/checkpoint_file.checkpoint.gz --stats \
## /path/to/stats_file.stats.gz --grid /path/to/grid_file.grid.gz --adapt_start 9 \
## -L 0.0001
## @endverbatim
##
## This will result in the following output:
##
## @verbatim
## builder = LearnerBuilder()
## builder.buildRegressor()
## builder.withTrainingDataFromARFFFile('/path/to/dataset_train.arff.gz')
## builder.withTestingDataFromARFFFile('/path/to/dataset_test.arff.gz')
## builder = builder.withGrid()
## builder.fromFile('/path/to/grid_file.grid.gz')
## builder.withStartingIterationNumber(9)
## builder = builder.withSpecification()
## builder.withAdaptPoints(400)
## builder.withLambda(0.00010000000000000000)
## builder.withIdentityOperator()
## builder = builder.withStopPolicy()
## builder.withGridSizeLimit(50000)
## builder = builder.withCGSolver()
## builder.withAccuracy(0.000100)
## builder.withImax(500)
## builder.withProgressPresenter(InfoToScreenRegressor())
## builder.withProgressPresenter(InfoToFile('/path/to/stats_file.stats.gz'))
## checkpointController = CheckpointController('checkpoint_file.checkpoint.gz', '/path/to')
## builder.withCheckpointController(checkpointController)
## @endverbatim
##
## You can then use this output to generate the Learner object.
# include all necessary classes
from pysgpp.extensions.datadriven.learner.LearnerBuilder import LearnerBuilder
from pysgpp.extensions.datadriven.controller import (CheckpointController,
                                                     InfoToScreenRegressor, InfoToFile)

## Paste the code as a string into your script:
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
checkpointController = CheckpointController('checkpoint_file.checkpoint.gz', '/path/to')
builder.withCheckpointController(checkpointController)"""

## Now you can execute the code and obtain the Learner object with
# create the Learner object
exec(code)
learner = builder.andGetResult()

## Usually if you want to restart your classifier from the checkpoint and to perform
## some more learning iterations, you had to repeat last iteration step even so you
## already have the alpha vector.
## Depending on the problem size this recomputation can be very expensive.
## Here is how you can omit it.
##
## Once the Learner object is loaded, e.g., as showed above,
## update results to recalculate errors per grid points and refine the grid
# update the values of errors per grid point
dataset = learner.dataContainer
trainSubset = dataset.getTrainDataset()
testSubset = dataset.getTestDataset()
learner.updateResults(learner.alpha, trainSubset, testSubset)
# refine the grid
learner.refineGrid()
# resume the learning process
learner.learnDataWithTest()

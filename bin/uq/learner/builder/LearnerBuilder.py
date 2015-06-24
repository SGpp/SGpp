##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)               #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################


from bin.data.ARFFAdapter import ARFFAdapter
from bin.data.CSVAdapter import CSVAdapter
from bin.data.DataContainer import DataContainer
from bin.learner.LearnedKnowledge import LearnedKnowledge
from bin.learner.TrainingSpecification import TrainingSpecification
from bin.learner.TrainingStopPolicy import TrainingStopPolicy
from bin.learner.Types import BorderTypes
from bin.learner.formatter import LearnedKnowledgeFormatter
from bin.learner.solver.CGSolver import CGSolver
from bin.uq.learner.Interpolant import Interpolant
from bin.uq.learner.Regressor import Regressor
from pysgpp import createOperationMultipleEval

from GridDescriptor import GridDescriptor
from RegressorSpecificationDescriptor import RegressorSpecificationDescriptor
import bin.utils.json as json


## Implement mechanisms to create customized learning system
#
# @section Examples Usage examples
#
# To create a learning system first define if it should be for classification
# @code
#import bin.learner.LearnerBuilder as LearnerBuilder
#builder = LearnerBuilder()
#builder = builder.buildClassifier()
# @endcode  
# or regression
# @code
#builder = builder.buildRegressor()
# @endcode
#
# LearnerBuilder is implementing <a href="http://en.wikipedia.org/wiki/Fluent_interface" target="parent">Fluent Interface design pattern</a>
# it means it operates as an automata, switching in some state
# where you can set all parameters associated with some category. For example to 
# define the grid parameters you switch the builder into GridDescriptor set with
# @code
#builder = builder.withGrid()...
# @endcode
# and then defines corresponding parameters:
# @code
#builder = builder.withGrid().withLevel(5).withBorder(Types.BorderTypes.TRAPEZOIDBOUNDARY)
# @endcode
# Builder can automatically switches to the next state
# @code
#builder.withGrid()...withCGSolver().withAccuracy(0.00000001)...
# @endcode
# After all parameters are set you can return the constructed learning system
# with
# @code
#builder.andGetResult()
# @endcode
# 
# The complete construction could look like following:
# @code
#classifier = builder.buildClassifier()\
#                     .withTrainingDataFromARFFFile("./datasets/classifier.train.arff")\
#                     .withTestingDataFromARFFFile("./datasets/classifier.test.arff")\
#                     .withGrid().withLevel(2)\
#                     .withSpecification().withLambda(0.00001).withAdaptPoints(2)\
#                     .withStopPolicy().withAdaptiveIterationLimit(1)\
#                     .withCGSolver().withImax(500)\
#                     .withProgressPresenter(InfoToFile("./presentor.test"))\
#                     .andGetResult()
# @endcode
#
# @section Parameters Parameters and where I can set them?
# @li <code>level</code>: <i>Gridlevel</i> - @link GridDescriptor.withLevel() <code>withGrid().withLevel(2)</code> @endlink
# @li <code>dim</code>: <i>Griddimension</i> - Dimension is identified from the %data set (by calling @link withTrainingDataFromARFFFile() <code>.withTrainingDataFromARFFFile(...)</code> @endlink)
# @li <code>adaptive</code>: <i>Using an adaptive Grid with NUM of refines</i> - @link StopPolicyDescriptor.withAdaptiveIterationLimit() <code>.withStopPolicy().withAdaptiveIterationLimit(10)</code> @endlink
# @li <code>adapt_points</code>: <i>Number of points in one refinement iteration</i> - @link SpecificationDescriptor.withAdaptPoints() <code>.withSpecification().withAdaptPoints(100)</code> @endlink
# @li <code>adapt_rate</code>: <i>Percentage of points from all refinable points in one refinement iteration</i> - @link SpecificationDescriptor.withAdaptRate() <code>.withSpecification().withAdaptRate(0.05)</code> @endlink
# @li <code>adapt_start</code>: <i>The index of adapt step to begin with</i> - Is know handled by loading the learner with specified iteration level from  CheckpointController using @link bin.controller.CheckpointController.CheckpointController.loadAll() <code>checkpointController.loadAll(10)</code> @endlink
# @li <code>adapt_threshold</code>: @link bin.learner.TrainingSpecification.TrainingSpecification.setAdaptThreshold() <i>refinement threshold</i> @endlink - @link SpecificationDescriptor.withAdaptThreshold() .withSpecification().withAdaptThreshold(0.003) @endlink
# @li <code>mode</code>: <i>Specifies the action to do</i> - Call corresponding method, i.e. @link Learner.Learner.applyData() applyData@endlink, @link Learner.Learner.learnData() learnData@endlink, @link Learner.Learner.learnDataWithTest() learnDataWithTest@endlink, @link  Learner.Learner.learnDataWithFolding() learnDataWithFolding@endlink
# @li <code>zeh</code>: <i>Specifies the action to do</i> - @link SpecificationDescriptor.withIdentityOperator() .withSpecification().withIdentityOperator()@endlink or @link SpecificationDescriptor.withLaplaceOperator() .withSpecification().withLaplaceOperator()@endlink
# @li <code>foldlevel</code>: <i>specifies the number of sets generated</i> - Is set in @link LearnerBuilder.LearnerBuilder.FoldingDescriptor FoldingDescriptor@endlink: <code> builder.withSequentialFoldingPolicy().withLevel(level)</code>
# @li <code>onlyfoldnum</code>: <i>Run only fold I in n-fold cross-validation. Default: run all</i> - @link bin.controller.CheckpointController.CheckpointController.generateFoldValidationJob() checkpointController.generateFoldValidationJob()@endlink generates a set of independent learners and a job script to run as SGE job array. In this way all or individual jobs can be ran either with SGE jobs or in console.
# @li <code>lambda</code>: <i>Lambda</i> - @link SpecificationDescriptor.withLambda() .withSpecification().withLambda(0.00001)@endlink
# @li <code>imax</code>: <i>Max number of iterations</i> - @link CGSolverDescriptor.withImax() .withCGSolver().withImax(500)@endlink
# @li <code>accuracy</code>: <i>Specifies the accuracy of the CG-Iteration</i> - @link CGSolverDescriptor.withAccuracy() .withCGSolver().withAccuracy(0.0001)@endlink
# @li <code>max_accuracy</code>: <i>If the norm of the residuum falls below ACCURACY, stop the CG iterations</i> - @link CGSolverDescriptor.withThreshold() .withCGSolver().withThreshold(0.0000000001)@endlink
# @li <code>%data</code>: <i>Filename for the Datafile.</i> - @link LearnerBuilder.withTrainingDataFromARFFFile() .withTestingDataFromARFFFile("./datasets/classifier.test.arff")@endlink
# @li <code>test</code>: <i>File containing the testdata</i> - @link LearnerBuilder.withTestingDataFromARFFFile() .withTestingDataFromARFFFile("./datasets/classifier.test.arff")@endlink
# @li <code>alpha</code>: <i>Filename for a file containing an alpha-Vector</i> -  <code>%learner = builder.andGetResult()\n learner.knowledge = LearnedKnowledgeFileAdapter().load("./alphas.arff")</code>
# @li <code>outfile</code>: <i>Filename where the calculated alphas are stored</i> - <code>@link LearnerBuilder.withProgressPresenter() .withProgressPresenter@endlink(@link bin.controller.InfoToFile.InfoToFile InfoToFile@endlink("./presentor.test"))</code>
# @li <code>gnuplot</code>: <i>In 2D case, the generated can be stored in a gnuplot readable format</i> - Some Graphs can now be plotted with @link bin.controller.InfoToGraph.InfoToGraph InfoToGraph@endlink
# @li <code>resolution</code>: <i>Specifies the resolution of the gnuplotfile</i> - Not used, as <code>gnuplot</code> is not yet implemented
# @li <code>stats</code>: <i>In this file the statistics from the test are stored</i> - Can be implemented as subclass from @link bin.controller.LearnerEventController.LearnerEventController LearnerEventController@endlink
# @li <code>polynom</code>: <i>Sets the maximum degree for basis functions</i> - @link GridDescriptor.withPolynomialBase() .withGrid().withPolynomialBase(2)@endlink
# @li <code>border</code> <i>Enables special border base functions</i> - @link GridDescriptor.withBorder() .withGrid().withBorder(Types.BorderTypes.TRAPEZOIDBOUNDARY)@endlink
# @li <code>trapezoid-boundary</code> <i>Enables boundary functions that have a point on the boundary for every inner point (Trapezoid)</i> - @link GridDescriptor.withBorder() .withGrid().withBorder(Types.BorderTypes.TRAPEZOIDBOUNDARY)@endlink
# @li <code>complete-boundary</code> <i>Enables boundary functions that have more points on the boundary than inner points</i> - @link GridDescriptor.withBorder() .withGrid().withBorder(Types.BorderTypes.COMPLETEBOUNDARY)@endlink
# @li <code>verbose</code>: <i>Provides extra output</i> - Set the suitable @link bin.controller.LearnerEventController.LearnerEventController LearnerEventController@endlink implementation, i.e. <code>.withProgressPresenter(InfoToScreen())</code>
# @li <code>normfile</code>: <i>For all modes that read %data via stdin. Normalizes %data according to boundaries in FILE</i> - Can be implemented as subclass of <code>@link bin.data.ARFFAdapter.ARFFAdapter ARFFAdapter@endlink</code>
# @li <code>reuse</code>: <i>Reuse alpha-values for CG</i> - @link bin.learner.LearnerBuilder.LearnerBuilder.CGSolverDescriptor.withAlphaReusing() <code>.withCGSolver().withAlphaReusing()</code>@endlink
# @li <code>seed</code>: <i>Random seed used for initializing</i> Is set in @link LearnerBuilder.LearnerBuilder.FoldingDescriptor FoldingDescriptor@endlink: <code>builder.withRandomFoldingPolicy().withSeed(level)</code>
# @li <code>regression</code>: <i>Use regression approach</i> - <code>@link buildRegressor() builder.buildRegressor()@endlink</code>
# @li <code>checkpoint</code>: <i>Filename for checkpointing</i> - <code> %controller = @link bin.controller.CheckpointController.CheckpointController.__init__ CheckpointController("classification_job")@endlink @link withCheckpointController() builder.withCheckpointController(controller)@endlink</code>
# @li <code>grid</code>: <i>Filename for Grid-resume</i> - <code>@link GridDescriptor.fromFile() .withGrid.fromFile("gridfile.gz")@endlink</code>
# @li <code>epochs_limit</code>: <i>Number of refinement iterations (epochs), MSE of test %data have to increase, before refinement will stop</i> - <code>@link StopPolicyDescriptor.withEpochsLimit() .withStopPolicy().withEpochsLimit(20)@endlink</code>
# @li <code>mse_limit</code> <i>If MSE of test %data fall below this limit, refinement will stop</i> - <code>@link StopPolicyDescriptor.withMSELimit() .withStopPolicy().withMSELimit(0.0003)@endlink</code>, also <code>@link StopPolicyDescriptor.withAccuracyLimit() ..withStopPolicy().withAccuracyLimit(0.95)@endlink</code> for classification accuracy
# @li <code>grid_limit</code> <i>If the number of points on grid exceed grid_limit, refinement will stop</i> - <code>@link StopPolicyDescriptor.withGridSizeLimit() .withStopPolicy().withGridSizeLimit(40000)@endlink</code>
#
class LearnerBuilder(object):

    def __init__(self):
        """
        Constructor
        """
        # created @link bin.learner.Learner.Learner Learner @endlink object
        self._learner = None
        self._gridDescriptor = None
        self._specificationDescriptor = None
        self._stopPolicyDescriptor = None

        # @link bin.controller.CheckpointController.CheckpointController
        # CheckpointController @endlink if any used
        self._checkpointController = None
        self._gridDescriptor = None

    def getLearner(self):
        """
        Returns the object of learner subclass, that is currently beeing
        constructed
        """
        return self._learner

    def getCheckpointController(self):
        """
        Returns the checkpoint controller
        """
        return self._checkpointController

    def buildRegressor(self):
        """
        Start building regressor
        """
        self._learner = Regressor()
        self._specificationDescriptor = RegressorSpecificationDescriptor(self)
        return self._specificationDescriptor

    def buildInterpolant(self):
        self._learner = Interpolant()
        return self

    def withGrid(self):
        """
        Start description of the grid
        """
        self._gridDescriptor = GridDescriptor(self)
        return self._gridDescriptor

    def withStartingIterationNumber(self, iteration):
        """
        Set the starting iteration number
        @param iteration: integer starting iteration number
        """
        self._learner.setCurrentIterationNumber(iteration)
        return self

    def withTrainingDataFromARFFFile(self, filename, name="train"):
        """
        Signals to use data from ARFF file for training dataset
        @param filename: Filename where to read the data from
        @param name: Category name, default: "train"
        """
        dataContainer = ARFFAdapter(filename).loadData(name)
        if self._learner.dataContainer is not None:
            dataContainer = self._learner.dataContainer.combine(dataContainer)
            self._learner.setDataContainer(dataContainer)
        else:
            self._learner.setDataContainer(dataContainer)
        return self

    def withTestingDataFromARFFFile(self, filename):
        """
        Signals to use data from ARFF file for testing dataset
        @param filename: Filename where to read the data from
        @return: LearnerBuilder object itself
        @todo (khakhutv) implement test for the method
        """
        adapter = ARFFAdapter(filename)
        dataContainer = adapter.loadData(DataContainer.TEST_CATEGORY)
        if self._learner.dataContainer is not None:
            dataContainer = self._learner.dataContainer.combine(dataContainer)
            self._learner.setDataContainer(dataContainer)
        else:
            self._learner.setDataContainer(dataContainer)
        return self

    def withTrainingDataFromCSVFile(self, filename, name="train"):
        """
        Signals to use data from CSV file for training dataset
        @param filename: Filename where to read the data from
        @param name: Category name, default: "train"
        """
        dataContainer = CSVAdapter(filename).loadData(name)
        if self._learner.dataContainer is not None:
            dataContainer = self._learner.dataContainer.combine(dataContainer)
            self._learner.setDataContainer(dataContainer)
        else:
            self._learner.setDataContainer(dataContainer)
        return self

    def withTestingDataFromCSVFile(self, filename):
        """
        Signals to use data from CSV file for testing dataset
        @param filename: Filename where to read the data from
        @return: LearnerBuilder object itself
        @todo (khakhutv) implement test for the method
        """
        adapter = CSVAdapter(filename)
        dataContainer = adapter.loadData(DataContainer.TEST_CATEGORY)
        if self._learner.dataContainer is not None:
            dataContainer = self._learner.dataContainer.combine(dataContainer)
            self._learner.setDataContainer(dataContainer)
        else:
            self._learner.setDataContainer(dataContainer)
        return self

    def withInitialAlphaFromARFFFile(self, filename):
        """
        Signals to use initial data for alpha vector from ARFF file
        @param filename: Filename where to read the data from
        @todo (khakhutv) implement test for the method
        """
        alpha = LearnedKnowledgeFormatter().deserializeFromFile(filename)
        self._learner.alpha = alpha
        self._learner.knowledge.setMemento(alpha)
        return self

    def withCheckpointController(self, controller):
        """
        Attaches checkpoint controller to the learner
        @param controller: Checkpoint controller which implements
        LearnerEventController
        """
        self._checkpointController = controller
        self._learner.attachEventController(self._checkpointController)
        self._checkpointController.setGrid(self._learner.grid)
        self._checkpointController.setLearnedKnowledge(self._learner.knowledge)
        self._checkpointController.setLearner(self._learner)
        return self

    def withProgressPresenter(self, presentor):
        """
        Attaches progress presentor to the learner
        @param presentor: progress presentor which implements
        LearnerEventController
        """
        self._learner.attachEventController(presentor)
        if self._learner.solver is not None:
            self._learner.solver.attachEventController(presentor)
        return self

    def andGetResult(self):
        """
        Returns the builded learner (regressor or interpolant),
        should be called in the and of construction
        """
        # @todo (khakhutv) construction of default parameters should be done
        # here
        if self._gridDescriptor is None:
            raise AttributeError('No grid is specified')
        if self._specificationDescriptor is None:
            raise AttributeError('Learner is not specified')
        return self._learner

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


from bin.learner.LearnedKnowledge import LearnedKnowledge
from bin.learner.Classifier import Classifier
from bin.learner.TrainingStopPolicy import TrainingStopPolicy

from bin.learner.TrainingSpecification import TrainingSpecification
from bin.learner.solver.CGSolver import CGSolver
from bin.learner.Types import BorderTypes

from bin.pysgpp import *
from bin.data.ARFFAdapter import ARFFAdapter
from bin.data.CSVAdapter import CSVAdapter
from bin.data.DataContainer import DataContainer
from bin.learner.Regressor import Regressor

import bin.utils.json as json
from bin.learner.folding.SequentialFoldingPolicy import SequentialFoldingPolicy
from bin.learner.folding.RandomFoldingPolicy import RandomFoldingPolicy
from bin.learner.folding.StratifiedFoldingPolicy import StratifiedFoldingPolicy
from bin.learner.folding.FilesFoldingPolicy import FilesFoldingPolicy
from bin.learner.formatter import LearnedKnowledgeFormatter

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
#                     .withStopPolicy().withAdaptiveItarationLimit(1)\
#                     .withCGSolver().withImax(500)\
#                     .withProgressPresenter(InfoToFile("./presentor.test"))\
#                     .andGetResult()
# @endcode
#
# @section Parameters Parameters and where I can set them?
# @li <code>level</code>: <i>Gridlevel</i> - @link GridDescriptor.withLevel() <code>withGrid().withLevel(2)</code> @endlink
# @li <code>dim</code>: <i>Griddimension</i> - Dimension is identified from the %data set (by calling @link withTrainingDataFromARFFFile() <code>.withTrainingDataFromARFFFile(...)</code> @endlink)
# @li <code>adaptive</code>: <i>Using an adaptive Grid with NUM of refines</i> - @link StopPolicyDescriptor.withAdaptiveItarationLimit() <code>.withStopPolicy().withAdaptiveItarationLimit(10)</code> @endlink
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
    
    ##created @link bin.learner.Learner.Learner Learner @endlink object
    __learner = None                  
    
    ##@link bin.controller.CheckpointController.CheckpointController 
    # CheckpointController @endlink if any used
    __checkpointController = None     
    __gridDescriptor = None 
    __specificationDescriptor = None
    __stopPolicyDescriptor = None
    __solverDescriptor = None
    
        
    ##
    # Default constuctor
    ##  
    def __init__(self):
        self.__learner = None
        self.__gridDescriptor = None 
        self.__specificationDescriptor = None
        self.__stopPolicyDescriptor = None
    
    
    ## Returns the object of learner subclass, that is currently beeing constructed
    # @return the object of learner subclass, that is currently beeing constructed
    def getLearner(self):
        return self.__learner

    
    ## Returns the checkpoint controller
    # @return the checkpoint controller
    def getCheckpointController(self):
        return self.__checkpointController



    ## Start building Regressor
    # @return: LearnerBuilder itself
    # @todo (khakhutv) implement test for building regressor
    ##         
    def buildRegressor(self):
        self.__learner = Regressor()
        return self.__buildCommonLearner(self.__learner)
    
    
    ## Start building Classifier
    # @return: LearnerBuilder itself
    ##
    def buildClassifier(self,):
        self.__learner = Classifier()
        return self.__buildCommonLearner(self.__learner)
    
    
    def __buildCommonLearner(self, learner):
        learnedKnowledge = LearnedKnowledge()
        learner.setLearnedKnowledge(learnedKnowledge)
        #stopPolicy = TrainingStopPolicy()
        #learner.setStopPolicy(stopPolicy)
        return self
    
    
    ## Start description of specification parameters for learner
    # @return: SpecificationDescriptor
    ##
    def withSpecification(self):
        self.__specificationDescriptor = LearnerBuilder.SpecificationDescriptor(self)
        return self.__specificationDescriptor


    ## Start description of parameters of CG-Solver for learner
    # @return: CGSolverDescriptor
    ##
    def withCGSolver(self):
        self.__solverDescriptor = LearnerBuilder.CGSolverDescriptor(self)
        return self.__solverDescriptor


    ## Start description of parameters of CG-Solver for learner
    # @return: GridDescriptor
    ##
    def withGrid(self):
        self.__gridDescriptor = LearnerBuilder.GridDescriptor(self)
        return self.__gridDescriptor


    ## 
    #Set the starting iteration number ane return the builder object
    #
    # @param iteration: integer starting iteration number
    # @return: LeanreBuilder
    ##
    def withStartingIterationNumber(self, iteration):
        self.__learner.setCurrentIterationNumber(iteration)
        return self
    

    ## 
    # Start description of parameters of stop-policy for learner
    # @return: StopPolicyDescriptor
    ##
    def withStopPolicy(self):
        self.__stopPolicyDescriptor = LearnerBuilder.StopPolicyDescriptor(self)
        return self.__stopPolicyDescriptor


    ## 
    # Returns the builded learner (regressor or classifier), should be called in the and of construction
    #
    # @return: Learner (Classifier of Regressor)
    ##
    def andGetResult(self):
        # @todo (khakhutv) construction of default parameters should be done here
        if self.__gridDescriptor == None:
            self.__gridDescriptor = LearnerBuilder.GridDescriptor()
        if self.__specificationDescriptor == None:
            self.__specificationDescriptor == LearnerBuilder.SpecificationDescriptor(self)
        if self.__learner.specification.getBOperator() == None:
            self.__learner.specification.setBOperator(
            createOperationMultipleEval(self.__learner.grid, self.__learner.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))

        return self.__learner


    ## 
    # Signals to use N-fold cross validation with sequential folding rule
    #
    # @return: FoldingDescriptor
    ##
    def withSequentialFoldingPolicy(self):
        self.__foldingPolicyDescroptor = LearnerBuilder.FoldingDescriptor(self, LearnerBuilder.FoldingDescriptor.SEQUENTIAL)
        return self.__foldingPolicyDescroptor
    
    
    ## 
    # Signals to use N-fold cross validation with random folding rule
    #
    # @return: FoldingDescriptor
    ##
    def withRandomFoldingPolicy(self):
        self.__foldingPolicyDescroptor = LearnerBuilder.FoldingDescriptor(self, LearnerBuilder.FoldingDescriptor.RANDOM)
        return self.__foldingPolicyDescroptor
    
    
    ## 
    # Signals to use N-fold cross validation with stratified folding rule
    #
    # @return: FoldingDescriptor
    ##
    def withStratifiedFoldingPolicy(self):
        self.__foldingPolicyDescroptor = LearnerBuilder.FoldingDescriptor(self, LearnerBuilder.FoldingDescriptor.STRATIFIED)
        return self.__foldingPolicyDescroptor
    
    
    ## 
    # Signals to use N-fold cross validation from a set of files
    #
    # @return: FoldingDescriptor
    ##
    def withFilesFoldingPolicy(self):
        self.__foldingPolicyDescroptor = LearnerBuilder.FoldingDescriptor(self, LearnerBuilder.FoldingDescriptor.STRATIFIED)
        return self.__foldingPolicyDescroptor

    
    ## 
    # Signals to use data from ARFF file for training dataset
    #
    # @param filename: Filename where to read the data from
    # @param name: Category name, default: "train"
    # @return: LearnerBuilder
    ##
    def withTrainingDataFromARFFFile(self, filename, name="train"):
        dataContainer = ARFFAdapter(filename).loadData(name)
        if self.__learner.dataContainer != None:
            self.__learner.setDataContainer(self.__learner.dataContainer.combine(dataContainer))
        else:
            self.__learner.setDataContainer(dataContainer)        
        return self
    

    ## 
    # Signals to use data from ARFF file for testing dataset
    #
    # @param filename: Filename where to read the data from
    # @return: LearnerBuilder object itself
    # @todo (khakhutv) implement test for the method
    def withTestingDataFromARFFFile(self, filename):
        dataContainer = ARFFAdapter(filename).loadData(DataContainer.TEST_CATEGORY)
        if self.__learner.dataContainer != None:
            self.__learner.setDataContainer(self.__learner.dataContainer.combine(dataContainer))
        else:
            self.__learner.setDataContainer(dataContainer)
        return self
    
    
    ## 
    # Signals to use data from CSV file for training dataset
    #
    # @param filename: Filename where to read the data from
    # @param name: Category name, default: "train"
    # @return: LearnerBuilder
    ##
    def withTrainingDataFromCSVFile(self, filename, name="train"):
        dataContainer = CSVAdapter(filename).loadData(name)
        if self.__learner.dataContainer != None:
            self.__learner.setDataContainer(self.__learner.dataContainer.combine(dataContainer))
        else:
            self.__learner.setDataContainer(dataContainer)        
        return self
    

    ## 
    # Signals to use data from CSV file for testing dataset
    #
    # @param filename: Filename where to read the data from
    # @return: LearnerBuilder object itself
    # @todo (khakhutv) implement test for the method
    def withTestingDataFromCSVFile(self, filename):
        dataContainer = CSVAdapter(filename).loadData(DataContainer.TEST_CATEGORY)
        if self.__learner.dataContainer != None:
            self.__learner.setDataContainer(self.__learner.dataContainer.combine(dataContainer))
        else:
            self.__learner.setDataContainer(dataContainer)
        return self
    
    
    ## 
    # Signals to use initial data for alpha vector from ARFF file
    #
    # @param filename: Filename where to read the data from
    # @return: LearnerBuilder object itself
    # @todo (khakhutv) implement test for the method
    def withInitialAlphaFromARFFFile(self, filename):
        alpha = LearnedKnowledgeFormatter().deserializeFromFile(filename)
        self.__learner.alpha = alpha
        self.__learner.knowledge.setMemento(alpha)
        return self


    ##
    # Attaches checkpoint controller to the learner
    #
    # @param controller: Checkpoint controller which implements LearnerEventController
    # @return: LearnerBuilder
    def withCheckpointController(self, controller):
        self.__checkpointController = controller
        self.__learner.attachEventController(self.__checkpointController)
        self.__checkpointController.setGrid(self.__learner.grid)
        self.__checkpointController.setLearnedKnowledge(self.__learner.knowledge)
        self.__checkpointController.setLearner(self.__learner)
        return self


    ##
    # Attaches progress presentor to the learner
    #
    # @param presentor: progress presentor which implements LearnerEventController
    # @return: LearnerBuilder
    def withProgressPresenter(self, presentor):
        self.__learner.attachEventController(presentor)
        if self.__learner.solver != None:
            self.__learner.solver.attachEventController(presentor)
        return self
        
        
        
    ## 
    # Grid Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the grid
    class GridDescriptor:
        __builder = None
        __deg = None
        __level = None
        __file = None
        __border = None
        __dim = None
        
        
        ## Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        def __init__(self, builder):
            self.__builder = builder
            self.__dim = self.__builder.getLearner().dataContainer.getDim()
            self.__deg = None
            self.__level = None
            self.__file = None
            self.__border = None
            
        
        ## 
        # Overrides built-in method
        # if method called is not a object method of this Descriptor, most probably it's a method of
        # LearnerBuilder so it tries to call the method from our builder
        #
        # @param attr: String for method name
        # @return: Method calling in LearnerBuilder
        def __getattr__(self, attr):
            grid = None
            if self.__file != None:
                gridFormatter = GridFormatter()
                grid = gridFormatter.deserializeFromFile(self.__file)
                self.__builder.getLearner().setGrid(grid)
            else:
                if self.__dim == None or self.__level == None:
                    raise AttributeError, "Not all attributes assigned to create grid"                
                if self.__border != None: 
                    if self.__border == BorderTypes.TRAPEZOIDBOUNDARY:
                        grid = Grid.createLinearTrapezoidBoundaryGrid(self.__dim)            
                    elif self.__border == BorderTypes.COMPLETEBOUNDARY:
                        grid = Grid.createLinearBoundaryGrid(self.__dim)            
                    else:
                        if self.__deg > 1:
                            grid = Grid.createModPolyGrid(self.__dim, self.__deg)
                        else:
                            grid = Grid.createModLinearGrid(self.__dim)
                else: #no border points
                        if self.__deg > 1:
                            grid = Grid.createPolyGrid(self.__dim, self.__deg)
                        else:
                            grid = Grid.createLinearGrid(self.__dim)
                            
                generator = grid.createGridGenerator()
                generator.regular(self.__level)
            self.__builder.getLearner().setGrid(grid)
            return getattr(self.__builder, attr)
            
        
        ## 
        # Defines the level of the grid
        #
        # @param level: level as integer
        # @return: GridDescriptor itself
        def withLevel(self, level):
            self.__level = level
            return self
        
        
        ##
        # Defines the polynomial base of the grid
        #
        # @param deg: degree of polynomial base as integer
        # @return: GridDescriptor itself
        ##
        def withPolynomialBase(self, deg):
            self.__deg = deg
            return self
        
        
        ##
        # Defines the border type of the grid
        #
        # @param type: border type as defined in bin.learner.Types.BorderTypes
        # @return: GridDescriptor itself
        ##
        def withBorder(self, type):
            self.__border = type
            return self
        
        
        ##
        # Indicates that grid should be restored from file
        #
        # @param filename: String name of file the grid should be restored from
        # @return: GridDescriptor itself
        ##
        def fromFile(self, filename):
            self.__file = filename
            return self
        
        
    ##
    # TrainingStopPolicy Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the training stop policy
    ##    
    class StopPolicyDescriptor:
        __builder = None
        __policy = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.__builder = builder
            self.__policy = TrainingStopPolicy()
            
        
        ##
        # Overrides built-in method
        # if method called is not a object method of this Descriptor, most probably it's a method of
        # LearnerBuilder so it tries to call the method from our builder
        #
        # @param attr: String for method name
        # @return: Method calling in LearnerBuilder
        ##    
        def __getattr__(self, attr):
            # if method called is not a object method of this Descriptor, most probably it's a method of
            # LearnerBuilder so we store results of descriptor and try to call the method from our builder 
            #if attr not in dir(self):
            #if none parameters are set, only one iteration has to be made
            if self.__policy.getAdaptiveIterationLimit()==None and \
            self.__policy.getAccuracyLimit() == None and \
            self.__policy.getEpochsLimit() == None and \
            self.__policy.getGridSizeLimit() == None and \
            self.__policy.getMSELimit() == None:
                self.__policy.setAdaptiveIterationLimit(0)
            self.__builder.getLearner().setStopPolicy(self.__policy)
            return getattr(self.__builder, attr)
            
        
        ##
        # Defines the maximal number of refinement steps
        #limit
        # @param limit: integer for maximal number of refinement steps
        # @return: StopPolicyDescriptor itself
        ##    
        def withAdaptiveItarationLimit(self, limit):
            self.__policy.setAdaptiveIterationLimit(limit)
            return self
         
        
        ##
        # Defines the maximal number of epochs MSE of test data can constantly increase
        #
        # @param limit: integer for maximal number of epochs
        # @return: StopPolicyDescriptor itself
        ##  
        def withEpochsLimit(self, limit):
            self.__policy.setEpochsLimit(limit)
            return self
         
        
        ##Defines the MSE for test data, which have to be arrived
        #
        # @param limit: float for MSE
        # @return: StopPolicyDescriptor itself
        ## 
        def withMSELimit(self, limit):
            self.__policy.setMSELimit(limit)
            return self
         
        
        ##
        # Defines the maximal number of points on grid
        #
        # @param limit: integer for maximal number of points on grid
        # @return: StopPolicyDescriptor itself
        ##  
        def withGridSizeLimit(self, limit):
            self.__policy.setGridSizeLimit(limit)
            return self
        
        
        ##
        # Defines the accuracy for test data, which have to be arrived
        #
        # @param limit: float for accuracy
        # @return: StopPolicyDescriptor itself
        ## 
        def withAccuracyLimit(self, limit):
            self.__policy.setAccuracyLimit(limit)
            return self
        
        
    
    ##
    # TrainingSpecification Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the training specification
    ##     
    class SpecificationDescriptor:
        __builder = None
        __specification = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.__builder = builder
            self.__specification = TrainingSpecification()
            
        
        ##
        # Overrides built-in method
        # if method called is not a object method of this Descriptor, most probably it's a method of
        # LearnerBuilder so it tries to call the method from our builder
        #
        # @param attr: String for method name
        # @return: Method calling in LearnerBuilder
        ##      
        def __getattr__(self, attr):
            # if method called is not a object method of this Descriptor, most probably it's a method of
            # LearnerBuilder so we try to call the method from our builder 
            if self.__specification.getCOperator() == None: #use laplace operator default
                self.__specification.setCOperator(createOperationLaplace(self.__builder.getLearner().grid))
                self.__specification.setCOperatorType('laplace')
            self.__builder.getLearner().setSpecification(self.__specification)
            return getattr(self.__builder, attr)
        
        
        ##
        # Specifies regression parameter of the learner
        #
        # @param value: float for regression parameter
        # @return: SpecificationDescriptor itself
        ##
        def withLambda(self, value):
            self.__specification.setL(value)
            return self
        
        ##
        # Specifies refinement threshold
        #
        # @param value: float for refinement threshold
        # @return: SpecificationDescriptor itself
        ##
        def withAdaptThreshold(self, value):
            self.__specification.setAdaptThreshold(value)
            return self
        
        
        ##
        # Specifies number of points, which have to be refined in refinement step
        #
        # @param value: integer for number of points to refine
        # @return: SpecificationDescriptor itself
        ##
        def withAdaptPoints(self, value):
            self.__specification.setAdaptPoints(value)
            return self
        
        
        ## Specifies rate from total number of points on grid, which should be refined
        #
        # @param value: float for rate 
        # @return: SpecificationDescriptor itself
        ##
        def withAdaptRate(self, value):
            self.__specification.setAdaptRate(value)
            return self
        
        
        ## Specifies to use laplace operator
        #
        # @return: SpecificationDescriptor itself
        ##
        def withLaplaceOperator(self, ):
            self.__specification.setCOperator(createOperationLaplace(self.__builder.getLearner().grid))
            self.__specification.setCOperatorType('laplace')
            return self
        
        
        ## Specifies to use identity operator
        #
        # @return: SpecificationDescriptor itself
        ##
        def withIdentityOperator(self, ):
            self.__specification.setCOperator(createOperationIdentity(self.__builder.getLearner().grid))
            self.__specification.setCOperatorType('identity')
            return self
        
    
    ##
    # CGSolver Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the CG-Solver
    ##     
    class CGSolverDescriptor:
        __builder = None
        __solver = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.__builder = builder
            self.__solver = CGSolver()
        
        
        ##
        # Overrides built-in method
        # if method called is not a object method of this Descriptor, most probably it's a method of
        # LearnerBuilder so it tries to call the method from our builder
        #
        # @param attr: String for method name
        # @return: Method calling in LearnerBuilder
        ##   
        def __getattr__(self, attr):
            # if method called is not a object method of this Descriptor, most probably it's a method of
            # LearnerBuilder so we try to call the method from our builder 
            self.__builder.getLearner().setSolver(self.__solver)
            return getattr(self.__builder, attr)
        
        
        ##
        # Defines the accuracy of CG-Solver
        #
        # @param accuracy: float for accuracy
        # @return: CGSolverDescriptor itself
        ##
        def withAccuracy(self, accuracy):
            self.__solver.setEpsilon(accuracy)
            return self
            
        
        ##
        # Defines the maxinmal number of iterations in CG algotihms
        #
        # @param imax: integer for maximal number of iteration in CG
        # @return: CGSolverDescriptor itself
        ##    
        def withImax(self, imax):
            self.__solver.setImax(imax)
            return self
        
        
        ##
        # Defines the maximal accuracy.
        # If the norm of the residuum falls below this threshold, stop the CG iterations
        #
        # @param threshold: maximal accuracy
        # @return: CGSolverDescriptor itself
        ##
        def withThreshold(self, threshold):
            self.__solver.setThreshold(threshold)
            return self
        
        
        ## The reusage of previous alpha data in the CG iteration
        # @return: CGSolverDescriptor itself
        def withAlphaReusing(self,):
            self.__solver.setReuse(True)
            return self
        
    
    ##
    # Folding Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning the usage for N-fold cross-validation
    ##   
    class FoldingDescriptor:
        
        SEQUENTIAL = 100 ## Sequential folding policy
        RANDOM = 200 ## Random folding policy
        STRATIFIED = 300 ## Stratified folding policy
        FILES = 400 ## Files folding policy
        
        __builder = None
        __level = None
        __type = None
        __policy = None
        __seed = None
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        # @param type: Type of folding policy that should be build
        ##
        def __init__(self, builder, type):
            self.__builder = builder
            self.__type = type
        
        
        ##
        # Overrides built-in method
        # if method called is not a object method of this Descriptor, most probably it's a method of
        # LearnerBuilder so it tries to call the method from our builder
        #
        # @param attr: String for method name
        # @return: Method calling in LearnerBuilder
        ##   
        def __getattr__(self, attr):
            # if method called is not a object method of this Descriptor, most probably it's a method of
            # LearnerBuilder so we try to call the method from our builder 
            if self.__builder.getLearner().dataContainer != None:
                dataContainer = self.__builder.getLearner().dataContainer
            else:
                raise Exception("Data not defined. Trainign data has to be defined before the folding policy")
            
            if self.__level == None:
                raise Exception("Folding level has to be defined")
            
            if self.__type == FoldingDescriptor.SEQUENTIAL:
                self.__policy = SequentialFoldingPolicy(dataContainer, self.__level)
            elif self.__type == FoldingDescriptor.RANDOM:
                self.__policy = RandomFoldingPolicy(dataContainer, self.__level, self.__seed)
            elif self.__type == FoldingDescriptor.STRATIFIED:
                self.__policy = StratifiedFoldingPolicy(dataContainer, self.__level, self.__seed)
            elif self.__type == FoldingDescriptor.FILES:
                self.__policy = FilesFoldingPolicy(dataContainer, self.__level)
            else:
                raise Exception("Folding type is not defined or is unproper")
            
            self.__builder.getLearner().setFoldingPolicy(self.__policy)
            return getattr(self.__builder, attr)
        
        
        ##
        # Defines the folding level
        #
        # @param level: integer folding level
        # @return: FoldingDescriptor itself
        ##
        def withLevel(self, level):
            self.__level = level
            
        ##
        # Defines the seed for random folding policy
        #
        # @param seed: integer seed
        # @return: FoldingDescriptor itself
        ##
        def withSeed(self, seed):
            self.__seed = seed
            
        

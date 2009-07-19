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

## @package LearnerBuilder
# @ingroup bin.learner
# @brief Realization of Builder Pattern for Learner Subclasses
# @version $CURR$

from bin.learner.LearnedKnowledge import LearnedKnowledge
from bin.learner.Classifier import Classifier
from bin.learner.TrainingStopPolicy import TrainingStopPolicy

from bin.learner.TrainingSpecification import TrainingSpecification
from bin.learner.CGSolver import CGSolver
from bin.learner.GridFileAdapter import GridFileAdapter
from bin.learner import Types

from bin.pysgpp import *
from bin.data.ARFFAdapter import ARFFAdapter
from bin.data.DataContainer import DataContainer
from bin.learner.Regressor import Regressor


class LearnerBuilder(object):
    learner = None
    gridDescriptor = None 
    specificationDescriptor = None
    stopPolicyDescriptor = None
    
        
    ##
    # Default constuctor
    ##  
    def __init__(self):
        self.learner = None
        self.gridDescriptor = None 
        self.specificationDescriptor = None
        self.stopPolicyDescriptor = None


    ##
    # Start building Regressor
    # @return: LearnerBuilder itself
    ##         
    def buildRegressor(self):
        self.learner = Regressor()
        #TODO: check if there is no better solution to deal with adapter
        learnedKnowledge = LearnedKnowledge()
        self.learner.setLearnedKnowledge(learnedKnowledge)
        stopPolicy = TrainingStopPolicy()
        self.learner.setStopPolicy(stopPolicy)
        return self
    
    
    ##
    # Start building Classifier
    # @return: LearnerBuilder itself
    ##
    def buildClassifier(self,):
        self.learner = Classifier()
        #TODO: check if there is no better solution to deal with adapter
        learnedKnowledge = LearnedKnowledge()
        self.learner.setLearnedKnowledge(learnedKnowledge)
        stopPolicy = TrainingStopPolicy()
        self.learner.setStopPolicy(stopPolicy)
        return self

    
    ##
    # Start description of specification parameters for learner
    # @return: SpecificationDescriptor
    ##
    def withSpecification(self):
        self.specificationDescriptor = LearnerBuilder.SpecificationDescriptor(self)
        return self.specificationDescriptor


    ##
    # Start description of parameters of CG-Solver for learner
    # @return: CGSolverDescriptor
    ##
    def withCGSolver(self):
        self.solverDescriptor = LearnerBuilder.CGSolverDescriptor(self)
        return self.solverDescriptor


    ##
    # Start description of parameters of CG-Solver for learner
    # @return: GridDescriptor
    ##
    def withGrid(self):
        self.gridDescriptor = LearnerBuilder.GridDescriptor(self)
        return self.gridDescriptor


    ##
    # Start description of parameters of stop-policy for learner
    # @return: StopPolicyDescriptor
    ##
    def withStopPolicy(self):
        self.stopPolicyDescriptor = LearnerBuilder.StopPolicyDescriptor(self)
        return self.stopPolicyDescriptor


    ##
    # Returns the builded learner (regressor or classifier),
    # should be called in the and of construction
    #
    # @return: Learner (Classifier of Regressor)
    ##
    def andGetResult(self):
        #TODO: construction of default parameters should be done here
        return self.learner


    ##
    # Signals to use N-fold cross validation
    # with sequential folding rule
    #
    # @return: FoldingDescriptor
    ##
    def withSequentialFoldingPolicy(self):
        return self

    
    ##
    # Signals to use data from ARFF file for training dataset
    # @param filename: Filename where to read the data from
    # 
    # @return: LearnerBuilder
    ##
    def withTrainingDataFromARFFFile(self, filename):
        dataContainer = ARFFAdapter(filename).loadData()
        if self.learner.dataContainer != None:
            self.learner.setDataContainer(self.learner.dataContainer.combine(dataContainer))
        else:
            self.learner.setDataContainer(dataContainer)        
        return self
    

    ##
    # Signals to use data from ARFF file for testing dataset
    #
    # @param filename: Filename where to read the data from
    # @return: LearnerBuilder
    ##
    def withTestingDataFromARFFFile(self, filename):
        dataContainer = ARFFAdapter(filename).loadData(DataContainer.TEST_CATEGORY)
        if self.learner.dataContainer != None:
            self.learner.setDataContainer(self.learner.dataContainer.combine(dataContainer))
        else:
            self.learner.setDataContainer(dataContainer)
        return self


    ##
    # Attaches checkpoint controller to the learner
    #
    # @param controller: Checkpoint controller which implements LearnerEventController
    # @return: LearnerBuilder
    ##
    def withCheckpointController(self, controller):
        self.learner.attachEventController(controller)
        return self


    ##
    # Attaches progress presentor to the learner
    # @param controller: progress presentor which implements LearnerEventController
    #
    # @return: LearnerBuilder
    ##
    def withProgressPresentor(self, presentor):
        self.learner.attachEventController(presentor)
        if self.learner.solver != None:
            self.learner.solver.attachEventController(presentor)
        return self
        
        
        
    ##
    # Grid Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the grid
    ##
    class GridDescriptor:
        builder = None
        deg = None
        level = None
        file = None
        border = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.builder = builder
            self.dim = self.builder.learner.dataContainer.getDim()
            self.deg = None
            self.level = None
            self.file = None
            self.border = None
            
        
        ##
        # Overrides built-in method
        # if method called is not a object method of this Descriptor, most probably it's a method of
        # LearnerBuilder so it tries to call the method from our builder
        #
        # @param attr: String for method name
        # @return: Method calling in LearnerBuilder
        ##
        def __getattr__(self, attr):
            grid = None
            if self.file != None:
                adapter = GridFileAdapter()
                grid = adapter.load(self.file)
                self.builder.learner.setGrid(grid)
            else:
                if self.dim == None or self.level == None:
                    raise AttributeError, "Not all attributes assigned to create grid"                
                if self.border != None: 
                    if self.border == Types.BorderTypes.USCALEDBOUNDARY:
                        grid = Grid.createLinearBoundaryUScaledGrid(self.dim)            
                    elif self.border == Types.BorderTypes.COMPLETEBOUNDARY:
                        grid = Grid.createLinearBoundaryGrid(self.dim)            
                    else:
                        if self.deg > 1:
                            grid = Grid.createModPolyGrid(self.dim, self.deg)
                        else:
                            grid = Grid.createModLinearGrid(self.dim)
                else: #no border points
                        if self.deg > 1:
                            grid = Grid.createPolyGrid(self.dim, self.deg)
                        else:
                            grid = Grid.createLinearGrid(self.dim)
                            
                generator = grid.createGridGenerator()
                generator.regular(self.level)
            self.builder.learner.setGrid(grid)
            return getattr(self.builder, attr)
            
        
        ##
        # Defines the level of the grid
        #
        # @param level: level as integer
        # @return: GridDescriptor itself
        ##
        def withLevel(self, level):
            self.level = level
            return self
        
        
        ##
        # Defines the polynomial base of the grid
        #
        # @param deg: degree of polynomial base as integer
        # @return: GridDescriptor itself
        ##
        def withPolynomialBase(self, deg):
            self.deg = deg
            return self
        
        
        ##
        # Defines the border type of the grid
        #
        # @param type: border type as defined in bin.learner.Types.BorderTypes
        # @return: GridDescriptor itself
        ##
        def withBorder(self, type):
            self.border = type
            return self
        
        
        ##
        # Indicates that grid should be restored from file
        #
        # @param filename: String name of file the grid should be restored from
        # @return: GridDescriptor itself
        ##
        def fromFile(self, filename):
            self.file = filename
            return self
        
        
    ##
    # TrainingStopPolicy Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the training stop policy
    ##    
    class StopPolicyDescriptor:
        builder = None
        policy = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.builder = builder
            self.policy = TrainingStopPolicy()
            
        
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
            self.builder.learner.setStopPolicy(self.policy)
            return getattr(self.builder, attr)
            
        
        ##
        # Defines the maximal number of refinement steps
        #
        # @param limit: integer for maximal number of refinement steps
        # @return: StopPolicyDescriptor itself
        ##    
        def withAdaptiveItarationLimit(self, limit):
            self.policy.setAdaptiveIterationLimit(limit)
            return self
         
        
        ##
        # Defines the maximal number of epochs MSE of test data can constantly increase
        #
        # @param limit: integer for maximal number of epochs
        # @return: StopPolicyDescriptor itself
        ##  
        def withEpochsLimit(self, limit):
            self.policy.setEpochsLimit(limit)
            return self
         
        
        ##
        # Defines the MSE for test data, which have to be arrived
        #
        # @param limit: float for MSE
        # @return: StopPolicyDescriptor itself
        ## 
        def withMSELimit(self, limit):
            self.policy.setMSELimit(limit)
            return self
         
        
        ##
        # Defines the maximal number of points on grid
        #
        # @param limit: integer for maximal number of points on grid
        # @return: StopPolicyDescriptor itself
        ##  
        def withGridSizeLimit(self, limit):
            self.policy.setGridSizeLimit(limit)
            return self
        
        
        ##
        # Defines the accuracy for test data, which have to be arrived
        #
        # @param limit: float for accuracy
        # @return: StopPolicyDescriptor itself
        ## 
        def withAccuracyLimit(self, limit):
            self.policy.setAccuracyLimit(limit)
            return self
        
        
    
    ##
    # TrainingSpecification Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the training specification
    ##     
    class SpecificationDescriptor:
        builder = None
        specification = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.builder = builder
            self.specification = TrainingSpecification()
            
        
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
            if self.specification.getCOperator() == None: #use laplace operator default
                self.specification.setCOperator(self.builder.learner.grid.createOperationLaplace())
            self.builder.learner.setSpecification(self.specification)
            return getattr(self.builder, attr)
        
        
        ##
        # Specifies regression parameter of the learner
        #
        # @param value: float for regression parameter
        # @return: SpecificationDescriptor itself
        ##
        def withLambda(self, value):
            self.specification.setL(value)
            return self
        
        
        ##
        # Specifies number of points, which have to be refined in refinement step
        #
        # @param value: integer for number of points to refine
        # @return: SpecificationDescriptor itself
        ##
        def withAdaptPoints(self, value):
            self.specification.setAdaptPoints(value)
            return self
        
        
        ## Specifies rate from total number of points on grid, which should be refined
        #
        # @param value: float for rate 
        # @return: SpecificationDescriptor itself
        ##
        def withAdaptRate(self, value):
            self.specification.setAdaptRate(value)
            return self
        
        
        ## Specifies to use laplace operator
        #
        # @param value: float for rate 
        # @return: SpecificationDescriptor itself
        ##
        def withLaplaceOperator(self, ):
            self.specification.setCOperator(self.builder.learner.grid.createOperationLaplace())
            return self
        
        
        ## Specifies to use identity operator
        #
        # @param value: float for rate 
        # @return: SpecificationDescriptor itself
        ##
        def withIdentityOperator(self, ):
            self.specification.setCOperator(self.builder.learner.grid.createOperationIdentity())
            return self
        
    
    ##
    # CGSolver Descriptor helps to implement fluid interface patter on python
    # it encapsulates functionality concerning creation of the CG-Solver
    ##     
    class CGSolverDescriptor:
        builder = None
        solver = None
        
        
        ##
        # Constructor
        #
        # @param builder: LearnerBuilder which creates this Descriptor
        ##
        def __init__(self, builder):
            self.builder = builder
            self.solver = CGSolver()
        
        
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
            self.builder.learner.setSolver(self.solver)
            return getattr(self.builder, attr)
        
        
        ##
        # Defines the accuracy of CG-Solver
        #
        # @param accuracy: float for accuracy
        # @return: CGSolverDescriptor itself
        ##
        def withAccuracy(self, accuracy):
            self.solver.setAccuracy(accuracy)
            return self
            
        
        ##
        # Defines the accuracy for test data, which have to be arrived
        #
        # @param limit: integer for maximal number of iteration in CG
        # @return: CGSolverDescriptor itself
        ##    
        def withImax(self, imax):
            self.solver.setImax(imax)
            return self


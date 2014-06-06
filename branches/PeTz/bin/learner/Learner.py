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



from bin.pysgpp import *
from bin.learner.solver.CGSolver import CGSolver
from bin.learner.folding.FoldingPolicy import FoldingPolicy
import bin.utils.json as json
from bin.learner.TrainingStopPolicy import TrainingStopPolicy
from bin.learner.TrainingSpecification import TrainingSpecification
from bin.data.DataContainer import DataContainer
import types


## The class implements the routines common for classifier and regressor.
# However, since some methods are defined and note implemented, this class is
# abstract and cannot be instantiated.
#
# The class also implements the subject of <a href="http://en.wikipedia.org/wiki/Observer_pattern" target="new">the observer
# design pattern</a>.
#
# @section Observer design pattern 
# To customize the processing of progress information in SGPP the observer pattern 
# is used. The classes that want to be informed about events should implement LearnerEvenController
# and subscribe by the instance of Learner subclass
# with <code>attachEventController()</code>. After some event of LearnerEvents
# arise, the Learner subclass object calls method <code>handleLearningEvent()</code> by 
# all subscribers. As subscribers get a reference to the Learner object, 
# they can retrieve the attributes of the learner and process the information.
#
# <i>Roles</i>
# - Subject: Learner
# - Concrete subject: e.g. Classifier or Regressor
# - Observer: LernerEventController
# - Concrete Observer: e.g. InfoToScreen
#
# Observer can also want to retrieve the process information from LinearSolver. 
# See documentation of@link bin.learner.solver.LinearSolver.LinearSolver LinearSolver@endlink for more information.
class Learner(object):
    
    ##list of object listening to the learning events
    eventControllers = None 
    
    ##DataContainer object with training (and maybe test) data
    dataContainer = None    
    
    ##TrainingStopPolicy object associated with this Learner
    stopPolicy = None       
    
    ##TrainingSpecification object associated with this Learner
    specification = None    
    
    ##Grid of the Learner
    grid = None             
    
    ##LearnedKnowledge where alpha is stored
    knowledge = None        
    
    ##Implementation of folding policy if training with folding is used
    foldingPolicy = None    
    
    ##LinearSolver object associated with this Learner
    solver = None           
    
    ##DMSystemMatrix object associated with this Learner
    linearSystem = None     
    
    ##Number of current iterations
    iteration = None      
    
    ##list of train accuracy values measured in refinement iteration
    trainAccuracy = None   
    
    ##list of test accuracy values measured in refinement iteration
    testAccuracy = None       
    
    ##DataVector with current alpha vector
    alpha = None            
    
    ##List of average training accuracy data over all refinement iterations
    trainingOverall = None    
    
    ##List of average training accuracy data over all refinement iterations
    testingOverall = None    
    
    ##List of numbers of point on grid for different refinement iterations
    numberPoints = None       
    
    __SERIALIZABLE_ATTRIBUTES = ['eventControllers', 'dataContainer', 
                                 'stopPolicy', 'specification', 'grid', 
                                 'knowledge','foldingPolicy', 'solver']
    
    
    ## Constructor
    def __init__(self):
        self.eventControllers = []
        self.trainAccuracy = []
        self.trainingOverall = []
        self.testAccuracy = []
        self.testingOverall = []
        self.numberPoints = []
        self.iteration = 0
        self.eventControllers = []


    ## Learn data from training data set and use validation data set to prevent overfitting
    #
    # @param dataset: DataContainer object with data sets, default value None (initialized data set used)
    # @return: DataVector of alpha
    def learnDataWithTest(self, dataset = None):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STARTED)
        self.specification.setBOperator(createOperationMultipleEval(self.grid,
                  self.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        
        if dataset == None: dataset = self.dataContainer
        
        #learning step
        trainSubset = dataset.getTrainDataset()
        #testpoint = data.allPoint\points
        #testvalues = data.allValues\values
        testSubset = dataset.getTestDataset()

        while True: #repeat until policy says "stop"
            self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_STARTED)

            self.alpha = self.doLearningIteration(trainSubset)

            #calculate avg. error for training and test data and avg. for refine alpha
            self.updateResults(self.alpha, trainSubset, testSubset)
            
            self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE)
         
            self.iteration += 1
            
            if(self.stopPolicy.isTrainingComplete(self)): break
            
            #refine grid
            self.refineGrid()
       
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_COMPLETE)
        return self.alpha
    
    
    ## Refines Grid
    # the function is not implemented here
    def refineGrid(self):
        raise NotImplementedError

    ## Calculate the value of the function for given points
    #
    # @param points: DataVector set of points
    # @return: DataVector values 
    def applyData(self, points):
        self.notifyEventControllers(LearnerEvents.APPLICATION_STARTED)
        # if learner is restored from checkpoint, you need to create new B Operator
        if self.specification.getBOperator() == None:
            self.specification.setBOperator(self.grid.createOperationB())
        size = points.getNrows()
        dim = points.getNcols()
        values = DataVector(size)
        row = DataVector(dim)
        for i in xrange(size):
            points.getRow(i, row)
            values[i] = self.grid.eval(self.knowledge.getAlphas(), row)
        self.notifyEventControllers(LearnerEvents.APPLICATION_COMPLETE)
        return values
    

    ## Simple data learning
    #
    # @return: DataVector of alpha
    def learnData(self):
        self.notifyEventControllers(LearnerEvents.LEARNING_STARTED)
        self.specification.setBOperator(createOperationMultipleEval(self.grid,
                    self.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        
        while True: #repeat until policy says "stop"
            self.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)
            #learning step
            self.alpha = self.doLearningIteration(self.dataContainer)
            
            #calculate avg. error for training and test data and avg. for refine alpha
            self.updateResults(self.alpha, self.dataContainer)
            self.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            self.iteration += 1
            if(self.stopPolicy.isTrainingComplete(self)): break
            #refine grid
            self.refineGrid()
       
        self.notifyEventControllers(LearnerEvents.LEARNING_COMPLETE)
        return self.alpha


    ## Learn data with cross-fold validation
    #
    # @return: list of DataVector alpha in different folds
    def learnDataWithFolding(self,):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_FOLDING_STARTED)
        self.specification.setBOperator(createOperationMultipleEval(self.grid,
                  self.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
     
        alphas = []
        for dataset in self.foldingPolicy:
            alphas.append(self.learnDataWithTest(dataset))
            
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_FOLDING_COMPLETE)
        return alphas
    
    
    ## Perform one learning step
    #
    # @param set: DataContainer training data set
    # @return: DataVector alpha vector
    def doLearningIteration(self, set):
        #initialize values
        self.linearSystem = DMSystemMatrix(self.grid,
                                           set.getPoints(),
                                           self.specification.getCOperator(),
                                           self.specification.getL())
        size =  self.grid.getStorage().size() 
        # Reuse data from old alpha vector increasing its dimension
        if self.solver.getReuse() and self.alpha != None:
            alpha = DataVector(self.alpha)
            alpha.resize(size)
        # Use new alpha vector
        else:
            alpha = DataVector(size)
            alpha.setAll( 0.0 )
        b = DataVector(size)
        self.linearSystem.generateb(set.getValues(), b)
        #calculates alphas
        self.solver.solve(self.linearSystem, alpha, b, self.solver.getReuse(), 
                          False, self.solver.getThreshold())
        return alpha


    ## Evaluate  accuracy 
    # this method is not implemented!
    # @param dataContainer: DataContainer data set
    # @param alpha: DataVector alpha-vector
    def evalError(self, dataContainer, alpha):
        raise NotImplementedError
    
    
    ## Update different statistics about training progress
    # this method is not implemented!
    # @param alpha: DataVector alpha-vector
    # @param trainSubset: DataContainer with training data
    # @param testSubset: DataContainer with validation data, default value: None
    def updateResults(self, alpha, trainSubset, testSubset = None):
        raise NotImplementedError
    
    
    ## Returns the number of current iteration
    #
    # @return: integer iteration number
    def getCurrentIterationNumber(self,):
        return self.iteration
    
    
    ## Sets the number of current iteration
    #
    # @param value: integer new iteration number
    def setCurrentIterationNumber(self, value):
        self.iteration = value
    
    
    ## Add observer to the list
    #
    # @param observer: LearnerEventController object
    def attachEventController(self, observer):
        if observer not in self.eventControllers: self.eventControllers.append(observer)

    
    ## Remove observer from the list
    #
    # @param observer: LearnerEventController object
    def detachEventController(self,observer):
        if observer in self.eventControllers: self.eventControllers.remove(observer)
    
    
    ## Notify all observers about the new event
    #
    # @param event: LearnerEvents event
    def notifyEventControllers(self, event):
        for controller in self.eventControllers:
            controller.handleLearningEvent(self, event)
    
    
    ## Setter for data container
    # @param container: the data container object
    # @return: leaner itself
    def setDataContainer(self, container):
        self.dataContainer = container
        return self

    
    ## Setter for grid
    # @param grid: the grid obejct
    # @return: leaner itself
    def setGrid(self, grid):
        self.grid = grid
        return self

    
    ## Setter for training specification
    # @param specification: the training specification object
    # @return: leaner itself
    def setSpecification(self, specification):
        self.specification = specification
        return self

    
    ## Setter for training stop policy
    # @param policy: the training stop policy object
    # @return: leaner itself
    def setStopPolicy(self, policy):
        self.stopPolicy = policy
        return self
    
    
    ## Setter for linear solver
    # @param solver: the linear solver object
    # @return: leaner itself
    def setSolver(self, solver):
        self.solver = solver
        return self


    ## Setter for learned knowledge
    # @param knowledge: the learned knowledge object
    # @return: leaner itself
    def setLearnedKnowledge(self, knowledge):
        self.knowledge = knowledge
        return self
    
    
    ## Setter for folding policy
    # @param policy: the folding policy object
    # @return: leaner itself
    def setFoldingPolicy(self, policy):
        self.foldingPolicy = policy
        return self
    
    
    ## Converts list of float values to string, where floats are written in exponential fromat
    #
    # @param list list of floats
    # @return the string representation of the list 
    def __listOfFloatsToString(self, list):
        result = '['
        for item in list[0:-1]:
            result += "%e"%item + ","
        result += "%e"%list[-1] + ']'
        return result
            
    
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object. 
    def toString(self):
        serializationString = "'module' : '" + self.__module__ + "',\n"
        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            
            #integers, dictionaries can serialized with str()
            if type(attrValue) in [types.IntType, types.DictType] and attrName.find("__") != 0: 
                serializationString += "'" + attrName + "'" + " : " + str(attrValue) + ",\n"
            
            #store floats in exponential format
            elif type(attrValue) == types.FloatType:
                serializationString += "'" + attrName + "'" + " : " + "%e"%attrValue + ",\n"
            
            #store list of floats in exponential format
            elif type(attrValue) == types.ListType:
                if len(attrValue)>0 and type(attrValue[0]) == types.FloatType:
                    serializationString += "'" + attrName + "'" + " : " + self.__listOfFloatsToString(attrValue) + ",\n"
                else:
                    serializationString += "'" + attrName + "'" + " : " + str(attrValue) + ",\n"
            
            # serialize strings with quotes    
            elif type(attrValue) == types.StringType and attrName.find("__") != 0:
                serializationString += "'" + attrName + "'" + " : '" + attrValue + "',\n"
                
            #serialize knowledge
            elif attrName in self.__SERIALIZABLE_ATTRIBUTES :
                # grid and knowledge are stored by checkpoint controller itself
                if attrName not in ['grid', 'knowledge'] and attrName != 'foldingPolicy':
#                    try: 
                    serializationString += "'" + attrName + "'" + " : " + attrValue.toString() + ",\n" 
#                    except Exception as detail:
#                        print "**********atrname**********:",attrName

        return "{" + serializationString.rstrip(",\n").replace("'",'"') + "}"
    
    
    ## Restores the attributes of a subclass of Learner from the json object with attributes.
    #
    # @param jsonObject A json object.
    # @return The restored subclass object of Learner.
    def fromJson(self, jsonObject):
        self.trainAccuracy = jsonObject['trainAccuracy']
        self.trainingOverall = jsonObject['trainingOverall']
        self.testAccuracy = jsonObject['testingOverall']
        self.numberPoints = jsonObject['numberPoints']
        self.iteration = jsonObject['iteration']
        self.dataContainer = DataContainer.fromJson(jsonObject['dataContainer'])
        self.stopPolicy = TrainingStopPolicy.fromJson(jsonObject['stopPolicy'])
        self.specification = TrainingSpecification.fromJson(jsonObject['specification'])
        self.solver = CGSolver.fromJson(jsonObject['solver'])
        return self
    
    
    ##Restores the state which is saved in the given memento
    #
    #@param memento the memento object
    def setMemento(self, memento):
        self.fromJson(memento)
        
    
    ##Creates a new memento to hold the current state
    #
    #@return a new memento
    def createMemento(self):
        # ok, it's weird now since I've wrote the toString() method before 
        # createMemento(). Correct would be to create an Json Object and then 
        # convert it to the string
        jsonString = self.toString()
        jsonObject = json.JsonReader().read(jsonString)
        return jsonObject


## Constants of different learning events
class LearnerEvents:
    ## Learning process is started
    LEARNING_STARTED = 100          
    
    ## Learning process is complete            
    LEARNING_COMPLETE = 200            
    
    ## Learning with k-fold cross-validation is started         
    LEARNING_WITH_FOLDING_STARTED = 300         
    
    ## Learning with k-fold cross-validation is complete
    LEARNING_WITH_FOLDING_COMPLETE = 400      
    
    ## Learning / refinement step is started  
    LEARNING_STEP_STARTED = 500                 
    
    ## Learning / refinement step is complete
    LEARNING_STEP_COMPLETE = 600                
    
    ## Learning with testing (validation) is started
    LEARNING_WITH_TESTING_STARTED = 700         
    
    ## Learning with testing (validation) is complete
    LEARNING_WITH_TESTING_COMPLETE = 800        
    
    ## Initial / refinement step of learning with testing (validation) is started
    LEARNING_WITH_TESTING_STEP_STARTED = 900    
    
    ## Initial / refinement step of learning with testing (validation) is complete
    LEARNING_WITH_TESTING_STEP_COMPLETE = 1000  
    
    ## Applying of grid on data is started
    APPLICATION_STARTED = 1100                  
    
    ## Applying of grid on data is complete
    APPLICATION_COMPLETE = 1200                 
    
    ## Refining of the grid
    REFINING_GRID = 1300                        
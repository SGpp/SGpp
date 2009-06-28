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

## @package Learner
# @ingroup learner
# @brief Abstract class for learning
# @version $CURR$

from bin.pysgpp import *
from CGSolver import CGSolver
from bin.learner.FoldingPolicy import FoldingPolicy

class Learner(object):
    eventControllers = []

    dataContainer = None
    accuracy = None
    stopPolicy = None
    specification = None
    grid = None
    knowledge = None
    foldingPolicy = None
    solver = None
    linearSystem = None
    
    iteration = 0
    
    trainAccuracy = []
    testAccuracy = []
    alpha = None #DataVector
    trainingOverall = []
    testingOverall = []
    numberPoints = []
    
    
    def __init__(self):
        self.trainAccuracy = []
        self.trainingOverall = []
        self.testAccuracy = []
        self.testingOverall = []
        self.numberPoints = []


    def learnDataWithTest(self, dataset = None):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STARTED)
        self.specification.setBOperator(self.grid.createOperationB())
        #self.specification.setCOperator(self.grid.createOperationLaplace())
        
        if dataset == None: dataset = self.dataContainer
        
        #learning step
        trainSubset = dataset.getTrainDataset()
        #testpoint = data.allPoint\points
        #testvalues = data.allValues\values
        testSubset = dataset.getTestDataset()

        while True: #repeat until policy says "stop"
            #@todo: better event name
            self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_STARTED)

            self.alpha = self.doLearningIteration(trainSubset)

            #calculate avg. error for training and test data and avg. for refine alpha
            self.updateResults(self.alpha, trainSubset, testSubset)
            
            self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE)
         
            self.iteration += 1
            
            if(self.stopPolicy.isTrainingComplete(self)): break
            
            #refine grid
            #@todo: here comes some average alpha
            self.refineGrid()
            
       
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_COMPLETE)
        return self.alpha
    
    def refineGrid(self):
        raise NotImplementedError

    def applyData(self, points):
        self.notifyEventControllers(LearnerEvents.APPLICATION_STARTED)
        size = points.getSize()
        dim = points.getDim()
        values = DataVector(size)
        row = DataVector(1, dim)
        for i in xrange(size):
            points.getRow(i, row)
            values[i] = self.grid.eval(self.alpha, row)
        self.notifyEventControllers(LearnerEvents.APPLICATION_COMPLETE)
        return values
    
#    def calcResult(self,):
#        pass

    def learnData(self):
        self.notifyEventControllers(LearnerEvents.LEARNING_STARTED)
        self.specification.setBOperator(self.grid.createOperationB())
        #@fixme: its not always a laplace operation
        #self.specification.setCOperator(self.grid.createOperationLaplace())
        
        while True: #repeat until policy says "stop"
            self.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)
            #learning step
            self.alpha = self.doLearningIteration(self.dataContainer)
            #eval Error for training data and append it to other in this iteration
            self.trainAccuracy.append(self.evalError(self.dataContainer, self.alpha))
            #calculate avg. error for training and test data and avg. for refine alpha
            self.updateResults(self.alpha, self.dataContainer)
            self.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            self.iteration += 1
            if(self.stopPolicy.isTrainingComplete(self)): break
            #refine grid
            self.refineGrid()
       
        self.notifyEventControllers(LearnerEvents.LEARNING_COMPLETE)
        return self.alpha


    #@todo: make it possible to run jobs parallel
    def learnDataWithFolding(self,):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_FOLDING_STARTED)
        self.specification.setBOperator(self.grid.createOperationB())
        #@fixme: its not always a laplace operation
        #self.specification.setCOperator(self.grid.createOperationLaplace())
     
        alphas = []
        #@todo: can be called concurrently
        for dataset in self.foldingPolicy:
            alphas.append(self.learnDataWithTest(dataset))
            
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_FOLDING_COMPLETE)
        return alphas
    
    
    def doLearningIteration(self, set):
        #initialize values
        self.linearSystem = LinearSystem(set.getPoints(),
                                     set.getValues(), 
                                     self.specification.getBOperator(), 
                                     self.specification.getCOperator(),
                                     self.grid,
                                     self.specification.getL())
        #calculates alphas
        alpha = self.solver.solve(self.linearSystem)
        return alpha


    def evalError(self, dataContainer, alpha):
        raise NotImplementedError
    
#    def addToRefine(self, alphas, refineAlphas, trainingSet):
#        raise NotImplementedError
    
    def updateResults(self, alpha, trainSubset, testSubset = None):
        raise NotImplementedError
            
    def getCurrentIterationNumber(self,):
        return self.iteration
    
    def attachEventController(self, observer):
        if observer not in self.eventControllers: self.eventControllers.append(observer)

    def detachEventController(self,observer):
        if observer in self.eventControllers: self.eventControllers.remove(observer)
        
    def notifyEventControllers(self, event):
        for controller in self.eventControllers:
            controller.handleLearningEvent(self, event)
    
    def setDataContainer(self, container):
        self.dataContainer = container
        return self

    def setGrid(self, grid):
        self.grid = grid
        return self

    def setSpecification(self, specification):
        self.specification = specification
        return self

    def setStopPolicy(self, policy):
        self.stopPolicy = policy
        return self

    def setSolver(self, solver):
        self.solver = solver
        return self

    def setLearnedKnowledge(self, knowledge):
        self.knowledge = knowledge
        return self
    
    def setFoldingPolicy(self, policy):
        self.foldingPolicy = policy
        return self
    
class LearnerEvents:
    LEARNING_STARTED = 100
    LEARNING_COMPLETE = 200
    LEARNING_WITH_FOLDING_STARTED = 300
    LEARNING_WITH_FOLDING_COMPLETE = 400
    LEARNING_STEP_STARTED = 500
    LEARNING_STEP_COMPLETE = 600
    LEARNING_WITH_TESTING_STARTED = 700
    LEARNING_WITH_TESTING_COMPLETE = 800
    LEARNING_WITH_TESTING_STEP_STARTED = 900
    LEARNING_WITH_TESTING_STEP_COMPLETE = 1000
    APPLICATION_STARTED = 1100
    APPLICATION_COMPLETE = 1200
    REFINING_GRID = 1300
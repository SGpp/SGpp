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

# @version $HEAD$
from pysgpp import createOperationLaplace, createOperationIdentity
from bin.learner import Classifier, Regressor, LearnerEvents, LearnedKnowledge
from bin.controller.LearnerEventController import LearnerEventController
#from bin.data.ARFFAdapter import ARFFAdapter

import bin.utils.json as json
from bin.learner.formatter import LearnedKnowledgeFormatter, GridFormatter, LearnerFormatter
from bin.controller.InfoToScreenRegressor import InfoToScreenRegressor
from bin.controller.InfoToScreen import InfoToScreen
from bin.controller.InfoToFile import InfoToFile
try:
    from bin.controller.InfoToGraph import InfoToGraph
except ImportError: pass
import gzip, copy


## Class for handling events for storing and restoring of checkpoints.
# Responsible for storing of checkpoints during learning process and restoring 
# of learning from some checkpoint.
#
# @section Design_Pattern Memento Design Pattern
# CheckpointController has an ability to save the states of Learner,
# its Grid and LearnedKnowledge for different iterations of training process.
# To insure flexibility and extensibility of the system the design pattern 
# <a href="http://en.wikipedia.org/wiki/Memento_pattern" target="new">Memento</a>
# with several modifications was used.
#
# <i>Roles</i>
# - Originator: Grid, LearnedKnowledge, Learner
# - Memento: Grid object itself for Grid, DataVector with alpha values for 
# LearnedKnowledge object, json structure for Learner obejects
# - Caretaker: CeckpointController
# 
# <i>Why should we use Memento pattern here?</i>
# - Simplification of inner structure of classes, since they don't need to organized
# and store different states by themselves;
# - Decoupling of objects of business logic and their representations allows to assure
# compatibility of between different versions and support of different input/output
# formats without modification of business classes;
# - Client or Caretaker don't have to know about the organization of originator to store/restore checkpoints.
# 
# <i>Enhancements</i>\n
# To be able able to use the states after the program has once terminated, we have
# to make memento objects persistent. For this purpose the Classes GridFormatter,
# LearnedKnowledgeFormatter and LearnerFormatter are used, which provides the
# functionality of serialization and desiralization of corresponding objects into
# file or any other stream, e.g. socket or terminal.
#
# For example, to save the Grid object in file, controller requests its memento object with 
# states and then calls method GridFormatter().serializeToFile(<Grid memento object>, <file name>).
# 
# @section Examples How to restore a checkpoint
# Supposed, you have already learned and stored checkpoints and now you have files
# checkpoint.6.arff.gz, checkpoint.6.grid.gz and checkpoint.6.learner.gz from the 6th 
# iteration stored in your working directory. Now you can
# - Load the grid and/or with:\n
# <code>grid = CheckpointController('checkpoint').laodGrid(6)\n
# knowledge = CheckpointController('checkpoint').laodLearnedKnowledge(6)</code>\n
# - or load the whole learner, i.e. to change learning parameters and continue the learning\n
# <code>learner =  CheckpointController('checkpoint').loadAll(6)</code>\n
# 
class CheckpointController(LearnerEventController):

    ## Checkpoint title
    title = None
    
    ## Path to the checkpoint files
    path = None
    
    ## Reference to the Grid object
    grid = None
    
    ## Reference to the LearnedKnowledge object
    knowledge = None
    
    ## Number of the last iteration
    lastIteration = None
    
    ## Interval, which determines how often checkpoints should be set
    interval = None
    
    ## Reference to the Learner object
    learner = None
    
    ## Fold number if cross-validation is used
    fold = None
    
    ##Constructor
    # save checkpoint files will have a name like title.iteration.[grid | learner | arff].gz
    #@param title: string title for checkpoints
    #@param path: string absolute path to the checkpoint files
    #@param interval: integer defines the number of iteration between saved checkpoints, e.g. if interval = 1, only states from 0,5th, 10th, 15th ... iteration will be saved
    #@param fold: the folding level, if n-fold cross-validation is used
    def __init__(self, title, path = None, interval = None, fold = None):
        self.title = title
        self.path = path if path != None else '.'
        self.interval = interval if interval != None else 1
        if fold != None:
            self.fold = fold
            
    
    ## Composes checkpoint file name from path title and iteration number
    #@param iteration integer iteration number
    #@param fold: the folding level, if n-fold cross-validation is used
    #@return string composed name
    def composeName(self, iteration = None, fold=None):
        result = ""
        if fold != None:
            result = self.path + "/" + self.title + "." + "fold_" + str(fold)
        elif self.fold != None:
            result = self.path + "/" + self.title + "." + "fold_" + str(self.fold)
        else:
            result = self.path + "/" + self.title 
            
        if iteration != None:
            result = result + "." + str(iteration)
        return result
    
    
    ## Saves current Grid, LearnedKnowledge and Learner objects
    # @param iteration Integer iteration number.
    def saveAll(self, iteration):
        self.saveGrid(iteration)
        self.saveLearnedKnowledge(iteration)
        self.saveLearner(iteration)
    
    
    ## Loads the Learner object with corresponding Grid and LearnedKnowledge
    #@param iteration Integer iteration number.
    #@return the Learner object.
    def loadAll(self, iteration):
        self.grid = self.loadGrid(iteration)
        self.knowledge = self.loadLearnedKnowledge(iteration)
        self.learner = self.__loadLearner(iteration)
        self.learner.setGrid(self.grid)
        self.learner.setLearnedKnowledge(self.knowledge)
        self.learner.attachEventController(self)
        # setting C operator here, since the knowledge about grid is needed
        # recreation of B operator is not necessary since it will be set in learnData()
        cOperatorType = self.__getCOperatorType(iteration)
        if cOperatorType == 'laplace':
            self.learner.specification.setCOperator(createOperationLaplace(self.learner.grid))
            self.learner.specification.setCOperatorType('laplace')

        elif cOperatorType == 'identity':
            self.learner.specification.setCOperator(createOperationIdentity(self.learner.grid))
            self.learner.specification.setCOperatorType('identity')

        else:
            raise Exception('C Operator type is unknown')
        
        for controller in self.__getLearnerEventControllers(iteration):
            if controller['module'] == 'bin.controller.InfoToScreenRegressor':
                self.learner.attachEventController(InfoToScreenRegressor())
            elif controller['module'] == 'bin.controller.InfoToScreen':
                self.learner.attachEventController(InfoToScreen())
            elif controller['module'] == 'bin.controller.InfoToFile':
                self.learner.attachEventController(InfoToFile(controller['filename']))
            elif controller['module'] == 'bin.controller.InfoToGraph':
                self.learner.attachEventController(InfoToGraph(controller['filename']))
                
        for controller in self.__getSolverEventControllers(iteration):
            if controller['module'] == 'bin.controller.InfoToScreenRegressor':
                self.learner.solver.attachEventController(InfoToScreenRegressor())
            elif controller['module'] == 'bin.controller.InfoToScreen':
                self.learner.solver.attachEventController(InfoToScreen())
            elif controller['module'] == 'bin.controller.InfoToFile':
                self.learner.solver.attachEventController(InfoToFile(controller['filename']))
            elif controller['module'] == 'bin.controller.InfoToGraph':
                self.learner.solver.attachEventController(InfoToGraph(controller['filename']))
                
                
        return self.learner

    
    ## Setter for the current Grid object
    # @param grid The Grid object.
    def setGrid(self, grid):
        self.grid = grid


    ## Setter for the current LearnedKnowledge object
    #@param knowledge: @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink object
    def setLearnedKnowledge(self, knowledge):
        self.knowledge = knowledge


    ## Loads sg.Grid from the checkpoint  of given iteration
    #@param iteration: integer iteration number
    #@return Grid object
    def loadGrid(self, iteration):
        gridFile = self.composeName(iteration) + ".grid.gz"
        grid = GridFormatter().deserializeFromFile(gridFile)
        return grid


    ## Loads knowledge from the checkpoint of given iteration
    #@param iteration: integer iteration number
    #@return @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink object
    def loadLearnedKnowledge(self, iteration):
        knowledgeFile = self.composeName(iteration) + ".arff.gz"
        knowledgeMemento = LearnedKnowledgeFormatter().deserializeFromFile(knowledgeFile)
        knowledge = LearnedKnowledge()
        knowledge.setMemento(knowledgeMemento)
        return knowledge


    ## Saves current Grid to the checkpoint with given iteration
    #@param iteration: integer iteration number
    def saveGrid(self, iteration):
        gridFilename = self.composeName(iteration) + ".grid.gz"
        # @todo (khakhutv) (low) is there a good reason for Grid to be a Memento object for itself?
        gridMemento = self.grid.createMemento()
        GridFormatter().serializeToFile(gridMemento, gridFilename)



    ## Saves current LearnedKnowldge object to the checkpoint with given iteration
    #@param iteration: integer iteration number
    def saveLearnedKnowledge(self, iteration):
        knowledgeFile = self.composeName(iteration) + ".arff.gz"
        knowledgeMemento = self.knowledge.createMemento()
        LearnedKnowledgeFormatter().serializeToFile(knowledgeMemento, knowledgeFile)
        
        
    ## Saves the current Learner object to the checkpoint with given iteration
    #@param iteration integer iteration number
    def saveLearner(self, iteration):
        learnerFilename = self.composeName(iteration) + ".learner.gz"
        learnerMemento = self.learner.createMemento()
        LearnerFormatter().serializeToFile(learnerMemento, learnerFilename)
        
        
    def __getCOperatorType(self, iteration):
        # read data from file
        learnerFilename = self.composeName(iteration) + ".learner.gz"
        learnerMemento = LearnerFormatter().deserializeFromFile(learnerFilename)
        
        # reconstruct the object
        cOperatorType = learnerMemento['specification']['_TrainingSpecification__cOperatorType'];
        return cOperatorType
    
    
    
    def __getLearnerEventControllers(self, iteration):
        # read data from file
        learnerFilename = self.composeName(iteration) + ".learner.gz"
        learnerMemento = LearnerFormatter().deserializeFromFile(learnerFilename)
        controllers = learnerMemento['eventControllers']
        return controllers
    
    
    
    def __getSolverEventControllers(self, iteration):
        # read data from file
        learnerFilename = self.composeName(iteration) + ".learner.gz"
        learnerMemento = LearnerFormatter().deserializeFromFile(learnerFilename)
        controllers = learnerMemento['solver']['eventControllers']
        return controllers
    
    
    ## Loads the Learner object to the checkpoint with given iteration
    # @param iteration integer iteration number
    #@return The Learner object.
    def __loadLearner(self, iteration):
        # read data from file
        learnerFilename = self.composeName(iteration, self.fold) + ".learner.gz"
        
        learnerMemento = LearnerFormatter().deserializeFromFile(learnerFilename)
        
        # reconstruct the object
        learnerType = learnerMemento['module']
        if "bin.learner.Classifier" in learnerType:
            learner = Classifier()
        elif "bin.learner.Regressor" in learnerType:
            learner = Regressor()
        else:
            raise Exception("module name unknown")
        learner.setMemento(learnerMemento)
        return learner
    
    
    ## Learning event @link LearnerEventController.handleLearningEvent handler routine @endlink of LearnerEventController
    def handleLearningEvent(self, subject, status):  
        if status == LearnerEvents.LEARNING_STEP_COMPLETE or status == LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE:
            if subject.iteration % self.interval == 0:
                self.setLearner(subject)
                self.saveAll(subject.iteration)
    
    
    ## Setter for current Learner object
    # @param learner The Learner object         
    def setLearner(self, learner):
        self.setGrid(learner.grid)
        self.setLearnedKnowledge(learner.knowledge)
        self.learner = learner
        
    
    ## Generates a SGE array job script for concurrent performing of cross-fold
    # validation computations. The script can be then lunched using
    # \<code\>qsub -t 1-XXX \<scriptname\>.sge.job
    # @param email String with email-address, the status information from SGE should be sent to
    # @todo (khakhutv) write a test for this method
    def generateFoldValidationJob(self, email=""):
        
        if not self.learner.iteration == None:
            iteration = self.learner.iteration
        else:
            iteration = None
        
        
        # for each data sub-set generate Learner objects,
        # so they can work with one data set
        # Saving of learner will automatically save the subset
        fold = 0
        currentFold = self.fold
        for dataset in self.learner.foldingPolicy:
            self.fold = fold
            # save grid and learned knowledge
            # every fold has its own grid and alpha vector
            if not self.grid == None:
                self.saveGrid(iteration)
            if not self.knowledge.getAlphas() == None:
                self.saveLearnedKnowledge(iteration)
            
            newLearner = copy.copy(self.learner)
            newLearner.setDataContainer(dataset)

            
            learnerFilename = self.composeName(iteration, fold) + ".learner.gz"
            learnerMemento = newLearner.createMemento()
            LearnerFormatter().serializeToFile(learnerMemento, learnerFilename)
            fold = fold + 1
        
        self.fold = currentFold
            
        # generate SGE script 
        # @todo (khakhutv) how to manage the case where job file and checkpoints are in different directories?
        script = """
#!/bin/bash
#$-N %s
#$-j y
#$-o %s
#$-M """ + email + """
#$-m e
#$-v SGPP,OPT_TMP,PATH,LD_LIBRARY_PATH
. /etc/profile
echo "from bin.controller.CheckpointController import CheckpointController\nlearner = CheckpointController('""" \
+ self.title + """', '.', 1, $SGE_TASK_ID).loadAll(0)\nlearner.learnDataWithTest()" | python

echo "JOB_ID: $JOB_ID"
echo "SGE_TASK_ID: $SGE_TASK_ID"
"""
        sgeJobFilename = self.composeName() + ".sge.job"
        sgeJobFile = open(sgeJobFilename, "w")
        sgeJobFile.write(script)
        sgeJobFile.close()

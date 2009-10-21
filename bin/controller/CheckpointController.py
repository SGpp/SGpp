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

from bin.learner import Classifier, Regressor, LearnerEvents, LearnedKnowledge
from bin.controller.LearnerEventController import LearnerEventController
from bin.data.ARFFAdapter import ARFFAdapter

import bin.utils.json as json
from bin.learner import LearnedKnowledgeFormatter, GridFormatter, LearnerFormatter
import gzip


## Class for handling events for storing and restoring of checkpoints.
# Responsible for storing of checkpoints during learning process and restoring 
# of learning from some checkpoint.
#
# @section Memento Design Pattern
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
class CheckpointController(LearnerEventController):

    title = ""
    path = ""
    grid = None
    knowledge = None
    lastIteration = None
    interval = None
    learner = None
    
    ##Constructor
    # save checkpoint files will have a name like title.iteration.[grid | learner | arff].gz
    #@param title: string title for checkpoints
    #@param path: string absolute path to the checkpoint files
    #@param interval: integer defines the number of iteration between saved checkpoints, e.g. if interval = 5, only states from 0,5th, 10th, 15th ... iteration will be saved
    def __init__(self, title, path = ".", interval = 1):
        self.title = title
        self.path = path
        self.interval = interval
            
    
    ## Composes checkpoint file name from path title and iteration number
    #@param iteration integer iteration number
    #@return string composed name
    def composeName(self, iteration):
        return self.path + "/" + self.title + "." + str(iteration)
    
    
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
    #@return @link sg.Grid Grid @endlink object
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
        # @todo (khakhutv) is there a good reason for Grid to be a Memento object for itself?
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

    
    ## Loads the Learner object to the checkpoint with given iteration
    # @param iteration integer iteration number
    #@return The Learner object.
    def __loadLearner(self, iteration):
        # read data from file
        learnerFilename = self.composeName(iteration) + ".learner.gz"
        
        learnerMemento = LearnerFormatter().deserializeFromFile(learnerFilename)
        
        # reconstruct the object
        learnerType = learnerMemento['module']
        if "bin.learner.Classifier" in learnerType:
            learner = Classifier().fromJson( learnerMemento )
        elif "bin.learner.Regressor" in learnerType:
            learner = Regressor().fromJson( learnerMemento )
        else:
            raise Exception("module name unknown")
        learner.setMemento(learnerMemento)
        
        return learner
    
    
    ## Learning event @link LearnerEventController.handleLearningEvent handler routine @endlink of LearnerEventController
    def handleLearningEvent(self, subject, status):  
        if status == LearnerEvents.LEARNING_STEP_COMPLETE or status == LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE:
            if subject.iteration % self.interval == 0:
                self.saveAll(subject.iteration)
    
    
    ## Setter for current Learner object
    # @param learner The Learner object         
    def setLearner(self, learner):
        self.setGrid(learner.grid)
        self.setLearnedKnowledge(learner.knowledge)
        self.learner = learner
           



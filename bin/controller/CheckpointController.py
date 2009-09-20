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

## @package CheckpointController
# @ingroup bin.controller
# @brief Class for handling events for storing and restoring of checkpoints
# @version $HEAD$

from bin.learner.Learner import LearnerEvents
from bin.controller.LearnerEventController import LearnerEventController
from bin.learner.GridFileAdapter import GridFileAdapter
from bin.data.ARFFAdapter import ARFFAdapter


## Responsible for storing of checkpoints during learning process and restoring of learning from some checkpoint
class CheckpointController(LearnerEventController):

    title = ""
    path = ""
    grid = None
    knowledge = None
    lastIteration = None
    myLearnedKnowledge = None
    gridAdapter = None
    knowledgeAdapter = None
    interval = None
    
    ##Constructor
    #@param title: string title for checkpoints
    #@param path: string absolute path to the checkpoint files
    #@param interval: integer defines the number of iteration between saved checkpoints
    def __init__(self, title, path = "./", interval = 1):
        self.title = title
        self.path = path
        self.interval = interval
        self.gridAdapter = GridFileAdapter()
            
    
    ## Composes checkpoint file name from path title and iteration number
    #@param iteration: integer iteration number
    #@return string composed name
    def composeName(self, iteration):
        return self.path + "/" + self.title + "." + str(iteration)
    
    
    ## Save all grid and learned knowledge
    def saveAll(self, iteration):
        self.saveGrid(iteration)
        self.saveLearnedKnowledge(iteration)

    
    ## Setter for grid attribute
    def setGrid(self, grid):
        self.grid = grid


    ## Setter for knowledge attribute
    #@param knowledge: @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink object
    def setLearnedKnowledge(self, knowledge):
        self.knowledge = knowledge


    ## Loads sg.Grid from the checkpoint  of given iteration
    #@param iteration: integer iteration number
    #@return @link sg.Grid Grid @endlink object
    def loadGrid(self, iteration):
        gridFile = self.composeName(iteration) + ".grid.gz"
        grid = self.gridAdapter.load(gridFile)
        return grid


    ## Loads knowledge from the checkpoint of given iteration
    #@param iteration: integer iteration number
    #@return @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink object
    def loadLearnedKnowledge(self, iteration):
        knowledgeFile = self.composeName(iteration) + ".arff.gz"
        adapter = ARFFAdapter(knowledgeFile)
        knowledge = adapter.loadData()
        return knowledge


    ## Save Grid to the checkpoint with given iteration
    #@param iteration: integer iteration number
    def saveGrid(self, iteration):
        gridFilename = self.composeName(iteration) + ".grid.gz"
        self.gridAdapter.save(self.grid, gridFilename)


    ## Save knowledge to the checkpoint with given iteration
    #@param iteration: integer iteration number
    def saveLearnedKnowledge(self, iteration):
        knowledgeFile = self.composeName(iteration) + ".arff.gz"
        self.knowledgeAdapter = ARFFAdapter(knowledgeFile)
        self.knowledgeAdapter.save(self.knowledge.getAlphas())
    
    
    ## Learning event @link LearnerEventController.handleLearningEvent handler routine @endlink of LearnerEventController
    def handleLearningEvent(self, subject, status):  
        if status == LearnerEvents.LEARNING_STEP_COMPLETE or status == LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE:
            if subject.iteration % self.interval == 0:
                self.saveAll(subject.iteration)
           



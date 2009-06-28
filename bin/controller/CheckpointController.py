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
# @version $CURR$

from bin.learner.Learner import LearnerEvents
from bin.controller.LearnerEventController import LearnerEventController


class CheckpointController(LearnerEventController):
    """ generated source for CheckpointController

    """
    name = ""
    grid = None
    knowledge = None
    lastIteration = None
    myLearnedKnowledge = None

    def saveAll(self, iteration):
        pass

    def setGrid(self, grid):
        pass

    def setLearnedKnowledge(self, knowledge):
        pass

    def loadGrid(self, iteration):
        return

    def loadLearnedKnowledge(self, iteration):
        return

    def saveGrid(self, iteration):
        pass

    def saveLearnedKnowledge(self, iteration):
        pass
    
    def handleLearningEvent(self, subject, status):
        if status == LearnerEvents.LEARNING_STARTED:
            print "Dimension is:", subject.dataContainer.getDim()
            print "Number of datasets is:", subject.dataContainer.getSize()
            
        elif status == LearnerEvents.LEARNING_COMPLETE:
            print "Learning complete"
        
        elif status == LearnerEvents.LEARNING_WITH_FOLDING_STARTED:
            pass
        
        elif status == LearnerEvents.LEARNING_WITH_FOLDING_COMPLETE:
            print "Learning complete"
        
        elif status == LearnerEvents.LEARNING_STEP_STARTED:
            print "Adaptive Step: ", subject.iteration
            
        elif status == LearnerEvents.LEARNING_STEP_COMPLETE:
            print "Number of points: ", subject.numberPoints[-1]
            print "Correct classified on training data: ",subject.trainAccuracy[-1]
            
        elif status == LearnerEvents.LEARNING_WITH_TESTING_STARTED:
            print "Dimension is:", subject.dataContainer.getDim()
            print "Number of datasets is:", subject.dataContainer.getSize()
        
        elif status == LearnerEvents.LEARNING_WITH_TESTING_COMPLETE:
            print "Learning complete"
        
        elif status == LearnerEvents.LEARNING_WITH_TESTING_STEP_STARTED:
            print "Adaptive Step:", subject.iteration
            
        elif status == LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE:
            print "Number of points: ", subject.numberPoints[-1]
            print "Correct classified on training data: ",subject.trainAccuracy[-1]
            print "Correct classified on testing data:  ",subject.testAccuracy[-1]
            
        elif status == LearnerEvents.APPLICATION_STARTED:
            pass
        
        elif status == LearnerEvents.APPLICATION_COMPLETE:
            pass
        
        elif status == LearnerEvents.REFINING_GRID:
            print "Refining Grid"



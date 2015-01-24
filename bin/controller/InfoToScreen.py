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


from bin.learner.solver.LinearSolver import LinearSolverEvents
from bin.learner.Learner import LearnerEvents
from bin.controller.LearnerEventController import LearnerEventController
from bin.controller.SolverEventController import SolverEventController


## The class processes the progress information from Learner and LinearSolver and
# shows it on the terminal screen.
class InfoToScreen(LearnerEventController, SolverEventController):
    
    ##
    #Handles events from Linear Solver 
    #
    #@param subject: Linear Solver object
    #@param status: Event Status of type LinearSolverEvents
    ##
    def handleSolvingEvent(self, subject, status):
        if status == LinearSolverEvents.STARTING:
            print "Solving started"
        elif status == LinearSolverEvents.CALC_STARTING:
            print "Starting norm of residuum: %g" % subject.getResiduum()
        elif status == LinearSolverEvents.ITERATION_COMPLETE:
            print "delta: %g" % subject.getResiduum()
        elif status == LinearSolverEvents.COMPLETE:
            print "Solving Complete"
            #print "Number of iterations: %d (max. %d)" % (subject.getNumberIterations(), subject.getImax())
            print "Number of iterations: %d" % (subject.getNumberIterations())
            print "Final norm of residuum: %g" % subject.getResiduum()
            
            
    ##
    #Handles events from Learner 
    #
    #@param subject: Learner object
    #@param status: Event Status of type LearnerEvents
    ##        
    def handleLearningEvent(self, subject, status):
        if status == LearnerEvents.LEARNING_STARTED:
            print "Dimension is:", subject.dataContainer.getDim()
            print "Number of data entries is:", subject.dataContainer.getSize()
            
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
            print "Number of data entries is:", subject.dataContainer.getSize()
        
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

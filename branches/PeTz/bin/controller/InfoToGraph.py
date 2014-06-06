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


from bin.learner.LinearSolver import LinearSolverEvents
from bin.learner.Learner import LearnerEvents
from bin.controller.LearnerEventController import LearnerEventController
from bin.controller.SolverEventController import SolverEventController
import matplotlib
matplotlib.use('Agg') #to prevent the bug with matplot unable to open display
from matplotlib import pyplot
from matplotlib.lines import Line2D


## This class processes the information about the current state of the learning 
# process and presents it in form of a graph.
# In order to use this class the <a href="http://matplotlib.sourceforge.net" target="new">matplotlib library</a> is required.
class InfoToGraph(LearnerEventController, SolverEventController):
    
    ##The constant figure ID of the plot with learner information
    LEARNER_FIGURE = 1 
    
    ##The constant figure ID of the plot with solver information
    SOLVER_FIGURE = 2
        
    __solverFigure = None
    __learnerFigure = None
    __solverFilename = None
    __learnerFilename = None
    __solverResiduum = None
    __learnerResiduum = None
    __learnerTrainResiduum = None
    __learnerTestResiduum = None
    
    
    ## Constuctor
    # @param learnerFilename: file name where learner progress information has to be stored (eps file)
    # @param solverFilename: file name where solver progress information has to be stored (eps file)
    def __init__(self, learnerFilename = "learner_progress.eps", solverFilename="solver_progress.eps"):
        self.__solverFigure = pyplot.figure(self.SOLVER_FIGURE)
        self.__solverFigure.add_subplot(111)
        pyplot.ioff()
        self.__learnerFigure = pyplot.figure(self.LEARNER_FIGURE)
        self.__learnerFigure.add_subplot(111)
        pyplot.ioff()
        self.__learnerFilename = learnerFilename
        self.__solverFilename = solverFilename
        self.__solverResiduum = 0.0
        self.__learnerTrainResiduum = 0.0
        self.__learnerTestResiduum = 0.0
        
        
    ##
    #Handles events from Linear Solver 
    #
    #@param subject: Linear Solver object
    #@param status: Event Status of type LinearSolverEvents
    ##
    def handleSolvingEvent(self, subject, status):
        pyplot.figure(self.SOLVER_FIGURE)
        if status == LinearSolverEvents.CALC_STARTING:
            pyplot.clf()
            self.__solverResiduum = 0.0
            
        elif status == LinearSolverEvents.ITERATION_COMPLETE:
            iteration = subject.getNumberIterations()
            residuum = subject.getResiduum()
            line = Line2D([iteration-1, iteration], [self.__solverResiduum, residuum])
            ax = pyplot.axes()
            ax.add_line(line)
            ax.autoscale_view()
            ax.set_xlabel('Iteration')
            ax.set_ylabel('Residuum')
            pyplot.draw()
            pyplot.savefig(self.__solverFilename, format="eps")
            self.__solverResiduum = residuum
            
       
            
            
    ##
    #Handles events from Learner 
    #
    #@param subject: Learner object
    #@param status: Event Status of type LearnerEvents
    ##        
    def handleLearningEvent(self, subject, status):
        pyplot.figure(self.LEARNER_FIGURE)
        if status == LearnerEvents.LEARNING_STARTED:
            pyplot.clf()
            self.__learnerTrainResiduum = 0.0
            self.__learnerGridSize = 0.0

        elif status == LearnerEvents.LEARNING_STEP_COMPLETE:
            gridSize = subject.numberPoints[-1]
            residuum = subject.trainAccuracy[-1]
            line = Line2D([self.__learnerGridSize, gridSize], [self.__learnerTrainResiduum, residuum])
            ax = pyplot.axes()
            ax.add_line(line)
            ax.autoscale_view()
            ax.set_ylabel('Accuracy')
            ax.set_xlabel('Grid size')
            pyplot.draw()
            pyplot.savefig(self.__learnerFilename, format="eps")
            self.__learnerTrainResiduum = residuum
            self.__learnerGridSize = gridSize
            
        elif status == LearnerEvents.LEARNING_WITH_TESTING_STARTED:
            pyplot.clf()
            self.__learnerTrainResiduum = 0.0
            self.__learnerTestResiduum = 0.0
            self.__learnerGridSize = 0.0
            
        elif status == LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE:
            gridSize = subject.numberPoints[-1]
            trainResiduum = subject.trainAccuracy[-1]
            testResiduum = subject.testAccuracy[-1]
            line1 = Line2D([self.__learnerGridSize, gridSize], [self.__learnerTrainResiduum, trainResiduum],c='b')
            line2 = Line2D([self.__learnerGridSize, gridSize], [self.__learnerTrainResiduum, testResiduum], c='r')
        
            ax = pyplot.axes()
            ax.add_line(line1)
            ax.add_line(line2)
            ax.autoscale_view()
            ax.legend((line1,line2),('Train','Test'))
            ax.set_ylabel('Accuracy')
            ax.set_xlabel('Grid size')
            pyplot.draw()
            pyplot.savefig(self.__learnerFilename, format="eps")
            self.__learnerTrainResiduum = trainResiduum
            self.__learnerTestResiduum = testResiduum
            self.__learnerGridSize = gridSize


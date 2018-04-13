# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


## Abstract class for solution of system of linear equations.
#
# The class defines the methods that have to be implemented in class of with
# routines for solution of linear equations in order to be used with SGPP.
#
# The class also implements the subject of <a href="Observer_pattern" target="new">the observer
# design pattern</a>.
#
# @section LinearSolver_Observer design pattern
# To customize the processing of progress information in SGPP the observer pattern
# is used. The classes that want to be informed about events should implement SolverEvenController
# and subscribe by the instance of LinearSolver subclass with
# <code>attachEventController()</code>. After some event of LinearSolverEvents
# arise, LinearSolver object calls method <code>handleSolvingEvent()</code> by
# all subscribers. As subscribers get a reference to the LinearSolver object,
# they can retrieve the attributes of the solver and process the information.
#
# <i>Roles</i>
# - Subject: LinearSolver
# - Concrete subject: e.g. CGSolver
# - Observer: SolverEventController
# - Concrete Observer: e.g. InfoToScreen
#
# Observer can also want to retrieve the process information from Learner. See documentation of
# @link python.learner.Learner.Learner Learner@endlink for more information.
class LinearSolver(object):

    ##list of object listening to the solver events
    eventControllers = None

    ##Constructor
    def __init__(self):
        self.eventControllers = []


    ## Solver linear system
    # this method has to be implemented in concrete class
    #
    # @param linearSystem: DMSystemMatrix object of Linear System to solve
    def solve(self, linearSystem):
        raise NotImplementedError


    ## Add observer to the list
    #
    # @param observer: ProgressInfoPresentor object
    def attachEventController(self, observer):
        if observer not in self.eventControllers: self.eventControllers.append(observer)


    ## Remove observer from the list
    #
    # @param observer: ProgressInfoPresentor object
    def detachEventController(self, observer):
        if observer in self.eventControllers: self.eventControllers.remove(observer)


    ## Notify all observers about the new event
    #
    # @param event: LinearSolverEvents event
    def notifyEventControllers(self, event):
        for controller in self.eventControllers:
            controller.handleSolvingEvent(self, event)



## Constants of different solving events
class LinearSolverEvents(object):

    ## Solving is started
    STARTING = 100

    ## Solving is complete
    COMPLETE = 200

    ## An iteration of the iterative solver is complete
    ITERATION_COMPLETE = 300

    ## Calculation is started
    CALC_STARTING = 400

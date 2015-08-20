# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


## Abstract class of Subscribers of LinearSolverEvents. The classes that wants to obtain
# the progress notifications from LinearSolver subclasses should implement this class. See @link
# python.learner.solver.LinearSolver.LinearSolver documentation of Learner@endlink for details.
class SolverEventController(object):

    ##
    #Handles events from LinearSolver
    #
    #@param subject: LinearSolver object
    #@param event: Event Status of type LinearSolverEvents
    ##
    def handleSolvingEvent(self, subject, event):
        raise NotImplementedError()

    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.
    def toString(self):
        return "{'module' : '" + self.__module__ + "'}\n"

    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.
    def __repr__(self):
        return '{' + self.toString().lstrip("{").rstrip("}\n") + '}'

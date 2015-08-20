# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## Abstract class of Subscribers of LearnerEvents. The classes that wants to obtain
# the progress notifications from Learner should implement this class. See @link
# python.learner.Learner.Learner documentation of Learner@endlink for details.
class LearnerEventController(object):

    ##
    #Handles events from Learner
    #
    #@param subject: Learner object
    #@param event: Event Status of type LearnerEvents
    ##
    def handleLearningEvent(self, subject, event):
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

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


from pysgpp.extensions.datadriven.controller.InfoToScreen import InfoToScreen
from pysgpp.extensions.datadriven.learner.Learner import LearnerEvents

# @package python.controller
## Prints some regression specific information together with information
# processed by @link python.controller.InfoToScreen.InfoToScreen InfoToScreen @endlink
class InfoToScreenRegressor(InfoToScreen):



    ##
    #Handles events from Learner
    #
    #@param subject: Learner object
    #@param status: Event Status of type LearnerEvents
    ##
    def handleLearningEvent(self, subject, status):
        if status == LearnerEvents.LEARNING_STEP_COMPLETE:
            print "Number of points: ", subject.numberPoints[-1]
            print "MSE on training data: ",subject.trainAccuracy[-1]
            print "L2-norm of error: ", subject.getL2NormError()
            print "Min error: ", subject.getMinError()
            print "Max error: ", subject.getMaxError()

        elif status == LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE:
            print "Number of points: ", subject.numberPoints[-1]
            print "MSE on training data: ",subject.trainAccuracy[-1]
            print "MSE on testing data:  ",subject.testAccuracy[-1]
            print "L2-norm of error on training data: : ", subject.getL2NormError()
            print "Min error on training data: : ", subject.getMinError()
            print "Max error on training data: : ", subject.getMaxError()
        else:
            super(self.__class__, self).handleLearningEvent(subject, status)

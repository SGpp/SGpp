# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

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



from bin.controller.InfoToScreen import InfoToScreen
from bin.learner.Learner import LearnerEvents


## Prints some regression specific information together with information 
# processed by @link bin.controller.InfoToScreen.InfoToScreen InfoToScreen @endlink
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
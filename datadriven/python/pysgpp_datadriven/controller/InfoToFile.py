# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##############################################################################
                                    #
#############################################################################

from pysgpp.extensions.datadriven.controller.InfoToScreen import InfoToScreen
import sys

## The class processes the progress information from Learner and LinearSolver and
# stores it into the file.
class InfoToFile(InfoToScreen):
    
    ## Filename, where output should be written
    filename = None 
    
    ##
    #Constructor
    #
    #@param filename: filename where output should be written as string 
    def __init__(self, filename):
        self.filename = filename
        
        
    ##
    #Handles events from Linear Solver 
    #
    #@param subject: Linear Solver object
    #@param status: Event Status of type LinearSolverEvents
    ##       
    def handleSolvingEvent(self, subject, status):
        file = open(self.filename, "a")
        tmpout = sys.stdout
        sys.stdout = file
        InfoToScreen.handleSolvingEvent(self, subject, status)
        sys.stdout = tmpout
        file.close()
        
        
    ##
    #Handles events from Learner 
    #
    #@param subject: Learner object
    #@param status: Event Status of type LearnerEvents
    ##        
    def handleLearningEvent(self, subject, status):
        file = open(self.filename, "a")
        tmpout = sys.stdout
        sys.stdout = file
        InfoToScreen.handleLearningEvent(self, subject, status)
        sys.stdout = tmpout
        file.close()
        
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object.  
    def toString(self):
        serializationString = InfoToScreen.toString(self).rstrip().lstrip("{").rstrip("}")
        serializationString += ", 'filename' : '" + self.filename + "'"
        serializationString = "{" + serializationString + "}\n"
        return serializationString

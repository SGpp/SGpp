# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##############################################################################
                                    #
#############################################################################

from pysgpp import DataVector

## Class keeps all information, which was learned during the learning process.
# Currently, only the alpha vector is stored, by in the future more information may come.
class LearnedKnowledge(object):

    
    __alphas = None         #DataVector with learned alpha
    
    
    ## Constructor
    def __init__(self):
        self.__alphas = DataVector(1) #a dummy DataVector object
        

    ##Restores the state which is saved in the given memento
    #
    #@param memento the memento object
    def setMemento(self, memento):
        self.update(memento)
    
    
    ##Creates a new memento to hold the current state
    #
    #@return a new memento
    def createMemento(self):
        return self.__alphas


    ##Returns current alpha vector
    #
    #@return: DataVector of current alpha
    def getAlphas(self):
        return self.__alphas


    ##Alters the current alpha vector
    #
    #@param alpha: new alpha vector
    #@return: LearnedKnowledge itself
    def update(self, alpha):
        self.__alphas = alpha
        return self




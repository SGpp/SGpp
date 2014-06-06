##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007-2009 Dirk Plueger (Dirk.Pflueger@in.tum.de)            #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
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

from bin.utils.GzipSerializer import GzipSerializer

from bin.learner.LearnedKnowledge import LearnedKnowledge
from bin.pysgpp import DataVector, DataMatrix

## Provides functionality for the runtime serialization of the LearnedKnowledge object
#
# This design intends to separate the binary object representation and its
# business logic from the text representation that can be saved into file.
# The class is a part of <a href="http://en.wikipedia.org/wiki/Memento_pattern" 
# target="new">Memento design pattern</a> described in details in @link 
# bin.controller.CheckpointController.CheckpointController CheckpointController
# @endlink.
# 
# Currently, the LearnerKnowledge memento object is DataVector object with 
# alpha values.
class LearnedKnowledgeFormatter(GzipSerializer):


    ##Deserializes the LearnedKnowledge memento object from the file.
    # The file may or may be not gzip compressed.
    #@param filename The name of file with serialized Grid.
    #@return The LearnedKnowledgeMemento object.
    def deserializeFromFile(self, filename):
        if type(filename) != type(''): raise AttributeError, "Filename as destination expected"
        serializationStream = self.gzOpen(filename, "r")
        try:
            knowledge = self.deserialize(serializationStream)
        finally:
            serializationStream.close()
        return knowledge
    
    
    ##Deserializes the LearnedKnowledgeMemento object from the stream
    #
    #@param serializationStream The stream to deserialize.
    #@return The Grid object.
    def deserialize(self, serializationStream):
        alphas = self.__readAlphaARFF(serializationStream)
        return alphas
        
        
    ##Serializes the LearnedKnowledge memento object to the file
    #
    #@param memento LearnedKnowledgeMemento
    #@param filename The name of file where the Grid object should be serialized to.
    def serializeToFile(self, memento, filename):
        stream = self.gzOpen(filename, "w")
        try:
            self.serialize(memento, stream)
        finally:
            stream.close()
    
    
    ##Serializes the LearnedKnowledge memento object to the stream
    #
    #@param memento The LearnedKnowledgeMemento object.
    #@param serializationStream Output stream where the LearnedKnowledgeMemento object should be serialized to.
    def serialize(self, memento, serializationStream):
        self.__writeAlphaARFF(serializationStream, memento)
        
    
    ##Returns a string that represents the LearnedKnowledge memento object.
    #
    # @param memento The LearnedKnowledgeMemento object.
    # @return A string that represents the LearnedKnowledge memento object.        
    def toString(self, memento):
        return memento.toString()
    
    
    ## Reads alpha vector from the input stream
    #
    # @param serializationStream The stream.
    # @return alpha DataVector
    def __readAlphaARFF(self, serializationStream):
        data = self.__readDataARFF(serializationStream)
        alpha = DataVector(len(data["data"][0]))
        for i in xrange(len(data["data"][0])):
            alpha[i] = data["data"][0][i]
        return alpha
    
    
    ## Writes alpha vector to output stream
    #
    # @param fout The output stream
    # @param alpha DataVector of alpha
    def __writeAlphaARFF(self, fout, alpha):
        if hasattr(fout, "name"):
            fout.write("@RELATION \"%s ALPHAFILE\"\n\n" % fout.name)
        fout.write("@ATTRIBUTE alpha NUMERIC\n")
        
        fout.write("\n@DATA\n")
        
        for i in xrange(len(alpha)):
            fout.write("%1.20f\n" % alpha[i])
        
        
    ## Reads from an ARFF file stream:
    # The data is stored in lists. There is a value list for every dimension of the data set. e.g. 
    # [[2, 3],[1, 1]] are the data points P_1(2,1) and P_2(3,1)
    #
    # @param serializationStream An ARFF file stream
    # @return returns a set of a array with the data (named data), a array with the classes (named classes)
    def __readDataARFF(self, serializationStream):
        data = []
        classes = []
        hasclass = False
    
        # get the different section of ARFF-File
        for line in serializationStream:
            sline = line.strip().lower()
            if sline.startswith("%") or len(sline) == 0:
                continue
    
            if sline.startswith("@data"):
                break
            
            if sline.startswith("@attribute"):
                value = sline.split()
                if value[1].startswith("class"):
                    hasclass = True
                else:
                    data.append([])
        
        #read in the data stored in the ARFF file
        for line in serializationStream:
            sline = line.strip()
            if sline.startswith("%") or len(sline) == 0:
                continue
    
            values = sline.split(",")
            if hasclass:
                classes.append(float(values[-1]))
                values = values[:-1]
            for i in xrange(len(values)):
                data[i].append(float(values[i]))
                
        # return
        return {"data":data, "classes":classes}
        


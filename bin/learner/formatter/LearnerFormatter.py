// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

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



from bin.utils.GzipSerializer import GzipSerializer
from bin.learner import Classifier, Regressor
import bin.utils.json as json

## Provides functionality for the runtime serialization of the Learner subclasses.
#
# This design intends to separate the binary object representation and its
# business logic from the text representation that can be saved into file.
# The class is a part of <a href="http://en.wikipedia.org/wiki/Memento_pattern" 
# target="new">Memento design pattern</a> described in details in @link 
# bin.controller.CheckpointController.CheckpointController CheckpointController
# @endlink.
# 
# Currently, the Learner memento object is a dictionary with attributes describing 
# the Learner object.
class LearnerFormatter(GzipSerializer):
    
    
    ##Deserializes the Learner memento object from the stream.
    #
    #@param serializationStream The stream to deserialize.
    #@return The Learner memento object.
    def deserialize(self, serializationStream):
        text = GzipSerializer.deserialize(self, serializationStream)
        reader = json.JsonReader()
        jsonObject = reader.read(text)
        return jsonObject
    
    
    ##Serializes the Learner memento to the stream
    #
    #@param memento The Learner memento structure.
    #@param serializationStream output stream where memento should be serialized to.
    def serialize(self, memento, serializationStream):
        text = json.JsonWriter().write(memento)
        serializationStream.write(text)


    ##Serializes the Learner memento to the file
    #
    #@param memento The Learner memento structure
    #@param filename The name of file where memento should be serialized to.
    def serializeToFile(self, memento, filename):
        fstream = self.gzOpen(filename, "w")
        try:
            self.serialize(memento, fstream)
        finally:
            fstream.close()
            
            
    ##Returns a string that represents the Learner memento object.
    #
    # @param memento The LearnerMemento object.
    # @return A string that represents the Learner memento object. 
    def toString(self, memento):
        return json.JsonWriter().write(memento) 
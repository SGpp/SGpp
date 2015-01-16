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


from pysgpp import *

from bin.utils.GzipSerializer import GzipSerializer

## Provides functionality for the runtime serialization of the Grid object
#
# This design intends to separate the binary object representation and its
# business logic from the text representation that can be saved into file.
# The class is a part of <a href="http://en.wikipedia.org/wiki/Memento_pattern" 
# target="new">Memento design pattern</a> described in details in @link 
# bin.controller.CheckpointController.CheckpointController CheckpointController
# @endlink.
#
# However strict separation of grid object and its representation is not 
# implemented,as Formatter becomes from Grid object not a GridMemento object, 
# but a string serialization of Grid itself. So Grid plays the role of 
# GridMemento. The complete decoupling is a subject for future work.
#
# @todo Bring serialization and deserialization of the Grid object in GridMemento object
class GridFormatter(GzipSerializer):

    
    ##Deserializes the Grid memento object from the stream
    #
    #@param serializationStream The stream to deserialize.
    #@return The Grid memento object.
    def deserialize(self, serializationStream):
        text = serializationStream.read()
        return Grid.setMemento(text)
    
    
    ##Deserializes the Grid object from the file.
    # The file may or may be not gzip compressed.
    #@param filename The name of file with serialized Grid.
    #@return The Grid memento object.
    def deserializeFromFile(self, filename):
        grid = None
        fstream = self.gzOpen(filename, "r")
        try:
            grid = self.deserialize(fstream)
        finally:
            fstream.close()
        return grid


    ##Serializes grid to the stream
    #
    #@param memento: the Grid memento object
    #@param streamserializationStream: output stream where grid should be serialized to
    def serialize(self, memento, streamserializationStream):
        text = self.toString(memento)
        streamserializationStream.write(text)


    ##Serializes grid to the file
    #
    #@param memento: the Grid memento object
    #@param filename The name of file where the Grid object should be serialized to.
    #@todo (khakhutv) make the same variable names across similar functions
    def serializeToFile(self, memento, filename):
        fstream = self.gzOpen(filename, "w")
        try:
            self.serialize(memento, fstream)
        finally:
            fstream.close()
    
    
    ##Returns a string that represents the Grid object.
    #
    # @param memento The Grid memento object.
    # @return A string that represents the Grid object.        
    def toString(self, memento):
        return memento.serialize()
    
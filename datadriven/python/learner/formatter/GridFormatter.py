# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp import *

from pysgpp.extensions.datadriven.utils.GzipSerializer import GzipSerializer


## Provides functionality for the runtime serialization of the Grid object
#
# This design intends to separate the binary object representation and its
# business logic from the text representation that can be saved into file.
# The class is a part of <a href="http://en.wikipedia.org/wiki/Memento_pattern"
# target="new">Memento design pattern</a> described in details in @link
# python.controller.CheckpointController.CheckpointController CheckpointController
# @endlink.
#
# However strict separation of grid object and its representation is not
# implemented,as Formatter becomes from Grid object not a GridMemento object,
# but a string serialization of Grid itself. So Grid plays the role of
# GridMemento. The complete decoupling is a subject for future work.
#
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

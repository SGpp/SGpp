# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##############################################################################
                                    #
#############################################################################

from pysgpp.base import DataVector
from bin.learner.formatter import GridFormatter
import matplotlib.pyplot as plt 
from numpy import zeros, sqrt, ceil, floor


class GridImageFormatter(GridFormatter):
    '''
    classdocs
    '''


    ##Serializes grid to the stream
    #
    #@param memento: the Grid memento object
    #@param streamserializationStream: output stream where grid should be serialized to
    def serialize(self, memento, streamserializationStream):
        raise NotImplemented("The current method is not implemented in the class")
        
        
    ##Deserializes the Grid memento object from the stream
    #
    #@param serializationStream The stream to deserialize.
    #@return The Grid memento object.
    def deserialize(self, serializationStream):
        raise NotImplemented("The current method is not implemented in the class")
    
    
    ##Deserializes the Grid object from the file.
    # The file may or may be not gzip compressed.
    #@param filename The name of file with serialized Grid.
    #@return The Grid memento object.
    def deserializeFromFile(self, filename):
        raise NotImplemented("The current method is not implemented in the class")


    ##Serializes grid to the file
    #
    #@param memento: the Grid memento object
    #@param filename The name of file where the Grid object should be serialized to.
    #@todo (khakhutv) make the same variable names across similar functions
    def serializeToFile(self, memento, filename):
        fstream = self.gzOpen(filename, "w")
        
        try:
            figure = plt.figure()
            grid = memento
            storage = grid.getStorage()
            coord_vector = DataVector(storage.dim())
            points = zeros([storage.size(), storage.dim()])
            for i in xrange(storage.size()):
                point = storage.get(i)
                point.getCoords(coord_vector)
                points[i] = [j for j in coord_vector.array()]
            num_of_sublots = storage.dim()*(storage.dim()-1)/2
            rows = int(ceil(sqrt(num_of_sublots)))
            cols = int(floor(sqrt(num_of_sublots)))
            i = 1
            
            for x1 in xrange(1,storage.dim()):
                for x2 in xrange(2,storage.dim()+1):
                     figure.add_subplot(rows*100 + cols*10 + i)
                     figure.add_subplot(rows, cols, i)
                     plt.xlabel('x%d'%x1, figure=figure)
                     plt.ylabel('x%d'%x2, figure=figure)
                     plt.scatter(points[:,x1-1], points[:,x2-1], figure=figure)

                     i +=1
            plt.savefig(fstream, figure=figure)
            plt.close(figure)
        finally:
            fstream.close()
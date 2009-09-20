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

## @package DataContainer
# @ingroup bin.data
# @brief Container for points and corresponding function values 
# @version $HEAD$

from DataSpecification import DataSpecification
from bin.pysgpp import DataVector
from DataAdapter import DataAdapter
from DataEntry import DataEntry



class DataContainer(object):
    
    #Constants category names
    TRAIN_CATEGORY = 'train'
    TEST_CATEGORY = 'test'
    
    points = {}             #Dictionary for points from different categories of data sets
    values = {}             #Dictionary for values from different categories of data sets
    specification = None    #Specification of attributes of default data set
    dim = None              #Dimension of the default data set
    size = None             #Size of the default data set
    tempPoint = None        #Used for manipulations with points DataVector
    tempValue = None        #Used for manipulations with values DataVector
    name = None             #Category name of the default data set
        

    ##Implementation of iterator method next()
    #
    # @return: the next element in the container
    def next(self):
        for row in xrange(0,self.size):
            self.points[self.name].getRow(row, self.tempPoint)
            self.values[self.name].getRow(row, self.tempValue)
            yield DataEntry(self.tempPoint, self.tempValue)
        return
    
    
    ##Implementation of iterator method __getitem__()
    #
    # @param item: integer index of the item in container
    # @return: entry with given index in the container
    def __getitem__(self,item):
        self.points[self.name].getRow(item, self.tempPoint)
        self.values[self.name].getRow(item, self.tempValue)
        return DataEntry(self.tempPoint, self.tempValue)
    
    
    ##Implementation of iterator method __iter__()
    # iterates through the container
    def __iter__(self):
        return self.next()
    
    
    ##Returns the data set which belongs to certain category
    #
    # @param category: String category name (train or test)
    # @return: DataContainer only with requested data set 
    def getDataSubsetByCategory(self, category):
        result = DataContainer(DataVector(self.points[category]), DataVector(self.values[category]), category)
        return result
    
    
    ##Creates DataContainer with entries from the given list
    #
    # @param indices: list of indices
    # @param name: String for category name of data set, default: "train"
    # @return: DataContainer with entries from the given list
    def getDataSubsetByIndexList(self, indeces, name="train"):
        size = len(indeces)
        subset_points = DataVector(size, self.dim)
        subset_values = DataVector(size)
        row = DataVector(self.dim)
        points = self.getPoints()
        values = self.getValues()
        
        i = 0
        for index in indeces:
            points.getRow(index, row)
            subset_points.setRow(i, row)
            subset_values[i] = values[index]
            i = i + 1
        return DataContainer(subset_points, subset_values, name)
    
    
    ##Creates DataContainer only with train data set
    #
    # @return: DataContainer only with train data set
    def getTrainDataset(self):
        return self.getDataSubsetByCategory(self.TRAIN_CATEGORY)
    
    
    ##Creates DataContainer only with test data set
    #
    # @return: DataContainer only with test data set
    def getTestDataset(self):
        return self.getDataSubsetByCategory(self.TEST_CATEGORY)
    

    def load(self):
        return
    

    def normalize(self):
        return
    

#    def __init__(self, adapter):
#    def __init__(self, size, dim, name = "train"):

    ## Constructor
    # possible parameter combinations:
    # DataContainer(adapter)
    # DataContainer(size, dim, [name])
    # DataContianer(points, values, [name]
    #
    # @param adapter: Object implementing DataContainer
    # @param size: Integer size of data set
    # @param dim: Integer dimension of data set
    # @param name: category name, default: "train"
    # @param points: DataVector with points
    # @param values: DataVector with values
    def __init__(self, *args):
        self.points={}
        self.values={}
        if isinstance(args[0], DataAdapter): #takes (adapter: DataAdapter)
            pass
        elif type(args[0]) == type(0): #takes (size: int, dim: int, name="train")
            #FIXME: here IndexError is possible
            try:
                if args[2] is None:
                    self.name = self.TRAIN_CATEGORY
                else:
                    self.name = args[2]
            except IndexError:
                self.name = self.TRAIN_CATEGORY
                
            self.size = args[0]
            self.dim = args[1]
            self.points[self.name] = DataVector(self.size, self.dim)

            self.values[self.name] = DataVector(self.size)
            
            self.specification = DataSpecification()
        elif isinstance(args[0], DataVector): #takes (points: DataVector, values: DataVector, name="train")
            try:
                if args[2] is None:
                    self.name = self.TRAIN_CATEGORY
                else:
                    self.name = args[2]
            except IndexError:
                self.name = self.TRAIN_CATEGORY

            self.points[self.name] = args[0]
            self.values[self.name] = args[1]
            self.size = self.points[self.name].getSize()
            self.dim = self.points[self.name].getDim()
            
        self.tempPoint = DataVector(1,self.dim)
        self.tempValue = DataVector(1,1)
       

    ##Merge to Data container into one, so it stores several data sets with different categories
    #
    # @param container: DataContainer that has to be combined with the called one
    # @return: new DataContainer with several data sets
    def combine(self, container):
        newContainer = DataContainer(self.getPoints(), self.getValues(), self.name)
        return newContainer.__setSubContainer(container.getPoints(), container.getValues(), container.getName())
        
    
    ##Adds points and values into dictionaries
    #
    # @param points: DataVector new points
    # @param values: DataVector new values
    # @param name: String category name under which points and values should be stored
    # @return: DataContainer itself
    def __setSubContainer(self, points, values, name):
        self.points[name] = points
        self.values[name] = values
        return self
    

    ##Create DataVector of given size and dimension with 0 for all entries
    #
    # @param size: Integer size of the DataVector
    # @param dim: Integer dimension of the DataVector
    # @return: new DataVector with 0 for all entries
    def createNullVector(self, size, dim):
        vector = DataVector(size, dim)
        vector.setAll(0)
        return vector
    
    
    ## Returns points stored in the data set with default name
    #
    # @return: DataVector of points
    def getPoints(self):
        return self.points[self.name]
    
    
    ## Returns values stored in the data set with default name
    #
    # @return: DataVector of values
    def getValues(self):
        return self.values[self.name]
    
    
    ## Return the default name of data set
    #
    # @return: String category name
    def getName(self):
        return self.name
    
    
    ## Returns dimension of the default data set
    #
    # @return: Integer dimension
    def getDim(self):
        if self.dim == None:
            self.dim = self.points[self.name].getDim()
        return self.dim
    
    
    ## Returns size of the default data set
    #
    # @return: Integer size 
    def getSize(self):
        if self.size == None:
            self.size = self.points[self.name].getSize()
        return self.size
    
    
    ## Return tuple of points and values
    #
    # @return: tuple (DataVector points, DataVector values)
    def getPointsValues(self):
        return (self.getPoints(), self.getValues())



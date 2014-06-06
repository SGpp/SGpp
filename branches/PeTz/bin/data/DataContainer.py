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


from DataSpecification import DataSpecification
from bin.pysgpp import DataVector, DataMatrix
from DataAdapter import DataAdapter
from DataEntry import DataEntry
import ARFFAdapter
import types


## A collection of data
# It can contain data sets for different categories, like "train"
# and "test" data. Implements some common operation on sets, like combining of two
# containers to one, as well as access to data levels, e.g. subset, points, values.
#
# The objects of DataContainer are iterable, so user can iterate through the points and
# values of the subset of default category defined in the attribute @link DataContainer.name name@endlink.
class DataContainer(object):
    
    ##Constants category names - training data
    TRAIN_CATEGORY = 'train'
    
    ##Constants category names - validation data
    TEST_CATEGORY = 'test'
    
    ## Dictionary for points from different categories of data sets
    points = None             
    
    ##Dictionary for values from different categories of data sets
    values = None             
    
    ## Specification of attributes of default data set
    specifications = None  
    
    ## Dimension of the default data set
    dim = None        
    
    ##Size of the default data set      
    size = None             
    
    ##Used for manipulations with points DataVector
    tempPoint = None        
    
    ##Used for manipulations with values DataVector
    tempValue = None        
    
    ##Category name of the default data set
    name = None             
        

    ##Implementation of iterator method next()
    #
    # @return: the next element in the container
    def next(self):
        for row in xrange(0,self.size):
            self.points[self.name].getRow(row, self.tempPoint)
            self.tempValue = self.values[self.name][row]
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
    # @exception if requested category name doesn't exist
    def getDataSubsetByCategory(self, category):
        if self.points.has_key(category) and self.values.has_key(category):
            result = DataContainer(DataMatrix(self.points[category]), DataVector(self.values[category]), category)
            return result
        else:
            raise Exception("Requested category name doesn't exist")
    
    
    ##Creates DataContainer with entries from the given list
    #
    # @param indices: list of indices
    # @param name: String for category name of data set, default: "train"
    # @return: DataContainer with entries from the given list
    def getDataSubsetByIndexList(self, indices, name="train"):
        size = len(indices)
        subset_points = DataMatrix(size, self.dim)
        subset_values = DataVector(size)
        row = DataVector(self.dim)
        points = self.getPoints()
        values = self.getValues()
        
        i = 0
        for index in indices:
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
    

#    def __init__(self, adapter):
#    def __init__(self, size, dim, name = "train"):

    ## Constructor
    # possible parameter combinations:
    # DataContainer(adapter)
    # DataContainer(size, dim, [name="train", filename=None])
    # DataContianer(points, values, [name="train", filename=None])
    #
    # param adapter: Object implementing DataContainer
    # param size: Integer size of data set
    # param dim: Integer dimension of data set
    # param name: category name, default: "train"
    # param points: DataVector with points
    # param values: DataVector with values
    def __init__(self, *args):
        self.points={}
        self.values={}
        self.specifications = {} 
        try:
            if isinstance(args[0], DataAdapter): #takes (adapter: DataAdapter)
                # @todo (khakhutv) implement test for this part
                adapter = args[0]
                container = adapter.loadData()
                self.points = container.points
                self.values = container.values
                self.dim = container.dim
                self.size = container.size
                self.specifications = container.specifications
                self.name = container.name
            else:
                if type(args[0]) == types.IntType: #takes (size: int, dim: int, name="train")
                    try:
                        if args[2] is None:
                            self.name = self.TRAIN_CATEGORY
                        else:
                            self.name = args[2]
                    except IndexError:
                        self.name = self.TRAIN_CATEGORY
                        
                    self.size = args[0]
                    self.dim = args[1]
                    self.points[self.name] = DataMatrix(self.size, self.dim)
        
                    self.values[self.name] = DataVector(self.size)
                    specification = DataSpecification()
                    specification.createNumericAttributes(self.dim)
                    self.specifications[self.name] = specification
                    
                elif isinstance(args[0], DataMatrix): #takes (points: DataVector, values: DataVector, name="train", filename=None)
                    try:
                        if args[2] is None:
                            self.name = self.TRAIN_CATEGORY
                        else:
                            self.name = args[2]
                    except IndexError:
                        self.name = self.TRAIN_CATEGORY
                    
                    self.points[self.name] = args[0]
                    self.values[self.name] = args[1]
                    self.size = self.points[self.name].getNrows()
                    self.dim = self.points[self.name].getNcols()
                    
                    specification = DataSpecification()
                    specification.createNumericAttributes(self.dim)
                    
                    # if data comes from a file, note it in the specification
                    try:
                        if not args[3] is None:
                            specification.setFilename(args[3])
                            specification.setSaved()
                    except IndexError:
                        pass
                    
                    self.specifications[self.name] = specification
          
                self.tempPoint = DataVector(self.dim)
                self.tempValue = DataVector(1)
            
        except IndexError:
            raise Exception('Wrong or no attributes in constructor')


    ## Merge to Data container into one, so it stores several data sets with different categories
    #
    # @param container: DataContainer that has to be combined with the called one
    # @return: new DataContainer with several data sets
    def combine(self, container):
        newContainer = DataContainer(self.getPoints(), self.getValues(), self.name)
        for k in self.points.keys():
            if k != self.name:
                newContainer = newContainer.__setSubContainer(self.points[k], self.values[k], self.specifications[k], k)
        return newContainer.__setSubContainer(container.getPoints(), container.getValues(), container.getSpecifiction(), container.getName())
    
    
    ## Merges several data containers to one.
    # Unlike combine(), this method actually merges the set of points and values
    # and not just puts them to the different categories
    # @param cls python keyword (no specification required)
    # @param containerList list of DataContainer's
    @classmethod
    def merge(cls, containerList):
        if len(containerList) == 0:
            return None
        
        # determine the total number of entries
        size = 0
        for container in containerList:
            size += len(container.getValues())
            
        dim = container.getPoints().getNcols()
        
        # Copy data to the new DataVector's entry by entry
        allPoints = DataMatrix(size, dim)
        allValues = DataVector(size) 
        tmpVector = DataVector(dim)
        i = 0
        for container in containerList:
            points = container.getPoints()
            values = container.getValues()
            for j in xrange(len(values)):
                points.getRow(j, tmpVector)
                allPoints.setRow(i, tmpVector)
                allValues[i] = values[j]
                i += 1
        
        # return new DataContainer
        return DataContainer(allPoints, allValues)
            
    
    ##Adds points and values into dictionaries
    #
    # @param points: DataVector new points
    # @param values: DataVector new values
    # @param name: String category name under which points and values should be stored
    # @param specification specification
    # @return: DataContainer itself
    def __setSubContainer(self, points, values, specification, name):
        self.points[name] = points
        self.values[name] = values
        self.specifications[name] = specification
        return self
    

    ##Create DataVector of given size and dimension with 0 for all entries
    #
    # @param size: Integer size of the DataVector
    # @param dim: Integer dimension of the DataVector
    # @return: new DataVector with 0 for all entries
    def createNullVector(self, size, dim):
        vector = DataMatrix(size, dim)
        vector.setAll(0)
        return vector
    
    
    ## Returns points stored in the data set with default name
    # @param category String category name of the requested data ("train" or "test")
    # @return: DataVector of points
    def getPoints(self, category = None):
        if category == None:
            category = self.name
        return self.points[category]
    
    
    ## Returns values stored in the data set with default name
    # @param category String category name of the requested data ("train" or "test")
    # @return: DataVector of values
    def getValues(self, category = None):
        if category == None:
            category = self.name
        return self.values[category]
    
    ## Return the data specification of the default category
    # @return: the DataSpecification object
    def getSpecifiction(self):
        return self.specifications[self.name]
    
    
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
            self.size = self.points[self.name].getNrows()
        return self.size
    
    
    ## Return tuple of points and values
    #
    # @return: tuple (DataVector points, DataVector values)
    def getPointsValues(self):
        return (self.getPoints(), self.getValues())
    
    
    ##Returns a string that represents the object.
    #
    # @return A string that represents the object. 
    def toString(self):
        # save the data as a file, if it's not saved yet
        for category, specification in self.specifications.items():
            if not self.specifications[category].isSaved():
                ARFFAdapter.ARFFAdapter(self.specifications[category].getFilename())\
                        .save(self.getPoints(category), self.getValues(category), 
                        specification.getAttributes())
                specification.setSaved()
                
        serializedString = "'module' : '" + self.__module__ + "',\n"
        for category in self.specifications.keys():
            serializedString += "'" + category + "' : " + self.specifications[category].toString() + ",\n"
        return "{" + serializedString.rstrip(",\n") + "}\n"
    
    
    ## Restores the DataContainer object from the json object with attributes.
    #
    # @param cls python keyword for class method (no specification needed)
    # @param jsonObject A json object.
    # @return The restored DataContainer object.
    @classmethod
    def fromJson(cls, jsonObject):
        # initiate with train data, because they are always there
        specification = jsonObject['train']
        resultContainer = ARFFAdapter.ARFFAdapter(specification['filename']).loadData('train')
        # load data for other categories
        for category, specification in jsonObject.items():
            if not ( category == 'module' or category == 'train') :
                # @todo (khakhutv) implement test for the branch
                container = ARFFAdapter.ARFFAdapter(specification['filename']).loadData(category)
                resultContainer = resultContainer.combine(container)

        return resultContainer


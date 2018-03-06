import numpy as np

from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import copyGrid
from pysgpp import DataVector, Grid

import pysgpp.extensions.datadriven.uq.jsonLib as ju
import pysgpp.extensions.datadriven.utils.json as json
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCKnowledgeFormatter import ASGCKnowledgeFormatter


class ASGCKnowledge(object):
    """
    The ASGC knowledge class
    """

    def __init__(self):
        """
        Constructor
        """
        # {iteration: {qoi: <Grid>}}
        self.__grids = {}
        # {iteration: {qoi: {dtype: {t: <DataVector>}}}}
        self.__alphas = {}
        self.__iteration = 0


    @classmethod
    def initWithStandardValues(cls, grid, alpha):
        ans = ASGCKnowledge()
        ans.update(grid, alpha, "_", 0, KnowledgeTypes.SIMPLE, 0)
        return ans

    def getAvailableQoI(self):
        """
        get available quantities of interest
        @return: list of strings identifying the quantities of interest
        """
        if len(self.__alphas) == 0:
            raise Exception('No knowledge available')
        iteration = self.__alphas.iterkeys().next()
        return self.__alphas[iteration].keys()

    def getAvailableTimeSteps(self):
        """
        get available time steps
        @return: sorted list of floats
        """
        if len(self.__alphas) == 0:
            raise Exception('No knowledge available')
        iteration = self.__alphas.iterkeys().next()
        qoi = self.__alphas[iteration].iterkeys().next()
        dtype = self.__alphas[iteration][qoi].iterkeys().next()
        ts = self.__alphas[self.__iteration][qoi][dtype].keys()
        return sorted(ts)

    def getAvailableKnowledgeTypes(self):
        """
        @return list of available KnowledgeTypes
        """
        if len(self.__alphas) == 0:
            raise Exception('No knowledge available')
        iteration = self.__alphas.iterkeys().next()
        qoi = self.__alphas[iteration].iterkeys().next()
        return self.__alphas[iteration][qoi].keys()

    def getAvailableIterations(self):
        """
        get available iterations
        @return: sorted list of integes
        """
        return self.__grids.keys()

    def getIteration(self):
        """
        get current iteration number
        """
        return self.__iteration

    def setIteration(self, iteration):
        """
        set current iteration number
        """
        self.__iteration = iteration

    def hasAlpha(self, iteration, qoi, t, dtype):
        """
        Check if there is a coefficient vector for the given
        configuration.
        @param iteration: int iteration number
        @param qoi: string quantity of interest
        @param t: float time step
        @param dtype: KnowledgeType
        """
        return iteration in self.__alphas and \
            qoi in self.__alphas[iteration] and \
            dtype in self.__alphas[iteration][qoi] and \
            t in self.__alphas[iteration][qoi][dtype]

    def hasGrid(self, iteration, qoi):
        """
        Check if there is a grid available for the given configuration
        @param iteration: int iteration number
        @param qoi: string quantity of interest
        """
        return iteration in self.__grids and \
            qoi in self.__grids[iteration]

    def getGrid(self, qoi='_', iteration=None):
        """
        Get the grid for the given configuration
        @param qoi: string quantity of interest
        @param iteration: int, iteration number
        """
        if iteration is None:
            iteration = self.__iteration

        if self.hasGrid(iteration, qoi):
            return self.__grids[iteration][qoi]
        else:
            raise AttributeError('no grid for (i=%i, qoi=%s)' % (iteration,
                                                                 qoi))

    def getAlpha(self, qoi='_', t=0, dtype=KnowledgeTypes.SIMPLE,
                 iteration=None):
        """
        Get the coefficient vector for the given configuration
        @param qoi: string quantity of interest
        @param t: float time step
        @param dtype: KnowledgeType
        @param iteration: int, iteration number
        """
        if iteration is None:
            iteration = self.__iteration

        if self.hasAlpha(iteration, qoi, t, dtype):
            return self.__alphas[iteration][qoi][dtype][t]
        else:
            raise AttributeError('no knowledge for (i=%i, t=%g, qoi=%s, dtype=%i)' % (iteration, t, qoi,
                                                                                      dtype))

    def getAlphasByQoI(self, qoi='_', dtype=KnowledgeTypes.SIMPLE,
                       iteration=None):
        """
        Get all coefficient vectors for the given quantity of interest
        @param qoi: string quantity of interest
        @param iteration: int, iteration number
        """
        if iteration is None:
            iteration = self.__iteration

        if qoi in self.__alphas[iteration] and \
                dtype in self.__alphas[iteration][qoi]:
            return self.__alphas[iteration][qoi][dtype]
        else:
            raise AttributeError('no knowledge for (i=%i, qoi=%s, dtype=%i)' % (iteration, qoi,
                                                                                KnowledgeTypes.toString(dtype)))

    def getSparseGridFunction(self, qoi='_', t=0, dtype=KnowledgeTypes.SIMPLE,
                              iteration=None):
        """
        Get the sparse grid function (grid, alpha) for the given setting
        @param qoi: string quantity of interest
        @param t: float time step
        @param dtype: KnowledgeType
        @param iteration: int, iteration number
        """
        if iteration is None:
            iteration = self.__iteration

        # check if there is the desired grid
        if self.hasGrid(iteration, qoi):
            grid = self.getGrid(qoi, iteration=iteration)
        else:
            raise AttributeError('the grid for (i=%i, qoi=%s) does not exist' % (iteration, qoi))

        # check if there is the desired surplus vector
        if self.hasAlpha(iteration, qoi, t, dtype):
            alpha = self.getAlpha(qoi, t, dtype, iteration=iteration)
        else:
            raise AttributeError('the surplus vector for (i=%i, qoi=%s, t=%g, dtype=%i) does not exist'
                                 % (iteration, qoi, t, dtype))
        return grid, alpha

    def getAlphas(self):
        return self.__alphas

    def setAlphas(self, alphas):
        self.__alphas = alphas

    def getGrids(self):
        return self.__grids

    def setGrids(self, grids):
        self.__grids = grids

    def update(self, grid, alpha, qoi, t, dtype, iteration):
        """
        Update the knowledge
        @param grid: Grid
        @param alpha: numpy array surplus vector
        @param qoi: string quantity of interest
        @param t: float time step
        @param dtype: KnowledgeType
        @param iteration: int iteration number
        """
        # build dictionary
        if iteration not in self.__alphas:
            self.__alphas[iteration] = {}
            self.__grids[iteration] = {}

        if qoi not in self.__alphas[iteration]:
            self.__alphas[iteration][qoi] = {}
            self.__grids[iteration][qoi] = {}

        if dtype not in self.__alphas[iteration][qoi]:
            self.__alphas[iteration][qoi][dtype] = {}

        if t not in self.__alphas[iteration][qoi][dtype]:
            self.__alphas[iteration][qoi][dtype][t] = {}

        # store knowledge
        self.__iteration = iteration
        self.__grids[iteration][qoi] = grid.clone()
        self.__alphas[iteration][qoi][dtype][t] = np.array(alpha)

    def clearAlphas(self):
        self.__alphas = {}

    # ----------------------------------------------------------------
    # ASGCKnowledge File Formatter
    # ----------------------------------------------------------------
    def setMemento(self, memento):
        """
        Restores the state which is saved in the given memento
        @param memento: the memento object
        """
        self.fromJson(memento)

    def createMemento(self):
        """
        Creates a new memento to hold the current state
        """
        jsonString = self.toJson()
        jsonObject = json.JsonReader().read(jsonString)
        return jsonObject

    def writeToFile(self, filename):
        """
        Write knowledge object to file
        """
        m = self.createMemento()
        ASGCKnowledgeFormatter().serializeToFile(m, filename)

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the ASGC object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored ASGC object
        """
        knowledge = ASGCKnowledge()

        # restore iteration
        key = '_ASGCKnowledge__iteration'
        if key in jsonObject:
            knowledge.setIteration(int(jsonObject[key]))

        # restore surpluses: {iteration: {qoi: {dtype: {t: <Grid>}}}}
        key = '_ASGCKnowledge__grids'
        if key in jsonObject:
            grids = {}
            for iteration, v1 in jsonObject[key].items():
                d1 = {}
                for qoi, gridString in v1.items():
                    # undo the hack that made it json compatible
                    gridString = gridString.replace('__', '\n')\
                                           .encode('utf8')
                    # deserialize ...
                    grid = Grid.unserialize(gridString)
                    # ... and store it
                    d1[qoi] = grid
                grids[int(iteration)] = d1
            knowledge.setGrids(grids)

        # restore surpluses: {iteration: {qoi: {dtype: {t: <list float>}}}}
        key = '_ASGCKnowledge__alphas'
        if key in jsonObject:
            alphas = {}
            for iteration, v1 in jsonObject[key].items():
                d1 = {}
                for qoi, v2 in v1.items():
                    d2 = {}
                    for dtype, v3 in v2.items():
                        d3 = {}
                        for t, alpha in v3.items():
                            d3[float(t)] = DataVector(alpha).array()
                        d2[int(dtype)] = d3
                    d1[qoi] = d2
                alphas[int(iteration)] = d1
            knowledge.setAlphas(alphas)

        return knowledge

    def toJson(self):
        """
        @return: a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        attrName = "_ASGCKnowledge__alphas"
        attrValue = self.__getattribute__(attrName)
        serializationString += ju.parseAttribute(attrValue, attrName)

        attrName = "_ASGCKnowledge__grids"
        attrValue = self.__getattribute__(attrName)
        serializationString += ju.parseAttribute(attrValue, attrName)

        attrName = "_ASGCKnowledge__iteration"
        attrValue = self.__getattribute__(attrName)
        serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        # print "j-------------------------------------------"
        # print "{" + s + "}"
        # print "j-------------------------------------------"

        return "{" + s + "}"

    def __str__(self):
        return "%i, alphas=%i, grids=%i" % (self.__iteration,
                                            len(self.__alphas),
                                            len(self.__grids))

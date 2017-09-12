from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.uq_setting.UQSettingAdapter import UQSettingAdapter
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer

from ASGCStatistics import ASGCStatistics
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import copyGrid
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCKnowledge import ASGCKnowledge

import numpy as np
import sys


class ASGCUQManager(object):
    
    def __init__(self):
        # interacting objects
        self.uqSetting = None
        self.knowledge = ASGCKnowledge()
        self.sampler = None
        self.learner = None

        self.stats = ASGCStatistics()

        # test set
        self.testSet = None
        self.learnWithTest = False

        self.dataContainer = {}

        # simulation based parameters
        self.__params = None
        self._qoi = '_'
        self.__knowledgeTypes = [KnowledgeTypes.SIMPLE]
        self.__timeStepsOfInterest = [0]

        self.verbose = True

    def hasMoreSamples(self):
        return self.sampler.hasMoreSamples()

    def runNextSamples(self):
        samples = self.sampler.nextSamples(self.knowledge, self._qoi,
                                           self.__timeStepsOfInterest)
        if len(samples) > 0:
            self.uqSetting.runSamples(samples)
        self.learnData()

    # ----------------------------------------------------------------
    def __prepareDataContainer(self, data, name):
        """
        Prepare data for learning
        @param data: dictionary loaded from UQSetting
        @return dictionary {dtype: {t: <DataContainer>}}
        """
        ans = {}
        U = self.__params.activeParams().getIndependentJointDistribution()
        for dtype in self.getKnowledgeTypes():
            ans[dtype] = {}
            dim = self.sampler.getGrid().getStorage().getDimension()

            # prepare data container depending on the given knowledge type
            tmp = KnowledgeTypes.transformData(data, U, dtype)

            # load data for all time steps
            for i, (t, values) in enumerate(tmp.items()):
                size = len(values)
                mydata = np.ndarray((size, dim))
                sol = np.ndarray(size)
                for i, (sample, res) in enumerate(values.items()):
                    mydata[i, :] = np.array(sample.getActiveUnit())
                    sol[i] = float(res)
                ans[dtype][t] = DataContainer(points=mydata, values=sol, name=name)
        return ans

    # @profile
    def updateDataContainer(self, updateTestData=False):
        """
        Sets the training dataContainerDict container given a UQSetting

        WARNING: This method has severe performance issues. It needs
        to be improved such that it loads just the last computed
        chunk of samples.
        """
        # load the results for the new samples
        if len(self.dataContainer) == 0:
            ps = None
        else:
            # load time steps and quantity of interest
            # lookup for new samples
            ps = []
            dataContainer = self.dataContainer.itervalues().next().itervalues().next()
            for sample in self.uqSetting.getSamplesStats().values():
                if sample.getActiveUnit() not in dataContainer:
                    ps.append(sample)
#             print dataContainer.getSizeTrain(), " = ", self.uqSetting.getSize(), "-", len(ps)
#             assert dataContainer.getSizeTrain() == self.uqSetting.getSize() - len(ps)

        resultsDict = self.uqSetting.getTimeDependentResults(self.__timeStepsOfInterest, self._qoi, ps)
        # prepare the results as dictionary
        dataContainerDict = self.__prepareDataContainer(resultsDict, 'train')
        # set the new dataContainerDict container
        for dtype, values in dataContainerDict.items():
            if dtype not in self.dataContainer:
                self.dataContainer[dtype] = {}

            for t, newDataContainer in values.items():
                if t not in self.dataContainer[dtype]:
                    self.dataContainer[dtype][t] = newDataContainer
                else:
                    if newDataContainer.getSize() > 0:
                        self.dataContainer[dtype][t] = \
                            self.dataContainer[dtype][t].combine(newDataContainer)

        # if there is a test setting given, combine the train and the
        # test dataContainerDict container.
        # the test set does not change over time so update it just at the beginning
        if updateTestData and self.testSet is not None:
            dataContainerDict = self.testSet.getTimeDependentResults(self.__timeStepsOfInterest, self._qoi)
            dataContainerDict = self.__prepareDataContainer(dataContainerDict, 'test')
            for dtype, values in dataContainerDict.items():
                for t, newDataContainer in values.items():
                    if newDataContainer.getSize() > 0:
                        self.dataContainer[dtype][t] = \
                            self.dataContainer[dtype][t].combine(newDataContainer)

    def learnData(self):
        """
        Learn the available data
        """
        # learn the data
        self.updateDataContainer(updateTestData=self.sampler.getCurrentIterationNumber() == 1)
        if self.learnWithTest:
            self.learnDataWithTest()
        else:
            self.learnDataWithoutTest()

    # ----------------------------------------------------------------
    def learnDataWithoutTest(self, *args, **kws):
        # learn data
        if self.verbose:
            print "learning (i=%i, gs=%i, type=%s)" % (self.sampler.getCurrentIterationNumber(),
                                                       self.sampler.getGrid().getSize(),
                                                       self.sampler.getGrid().getTypeAsString())
        self.learner.grid = self.sampler.getGrid()
        for dtype, values in self.dataContainer.items():
            knowledge = {}
            if self.verbose:
                print KnowledgeTypes.toString(dtype)
            # do the learning
            for t, dataContainer in values.items():
                if self.verbose:
                    print "t = %g, " % t,
                sys.stdout.flush()
                if dataContainer is not None:
                    # learn data, if there is any available
                    self.learner.dataContainer = dataContainer
                    self.learner.learnData()

                    # update the knowledge
                    self.knowledge.update(copyGrid(self.learner.grid),
                                          self.learner.alpha,
                                          self._qoi,
                                          t,
                                          dtype,
                                          self.sampler.getCurrentIterationNumber() - 1)

                    # update results
                    self.stats.updateResults(dtype, t, self.learner)

            if self.verbose:
                print

    def learnDataWithTest(self, dataset=None, *args, **kws):
        if self.verbose:
            print "learning with test (i=%i, gs=%i)" % (self.sampler.getCurrentIterationNumber(),
                                                        self.sampler.getGrid().getSize())
        # learn data
        self.learner.grid = self.sampler.getGrid()
        for dtype, values in self.dataContainer.items():
            knowledge = {}
            # do the learning
            for t in np.sort(values.keys()):
                dataContainer = values[t]
                if self.verbose:
                    print "t = %g, " % t,
                sys.stdout.flush()
                if dataContainer is not None:
                    # learn data, if there is any available
                    self.learner.dataContainer = dataContainer
                    self.learner.learnDataWithTest(dtype=dtype)

                    # update the knowledge
                    self.knowledge.update(copyGrid(self.learner.grid),
                                          self.learner.alpha,
                                          self._qoi,
                                          t,
                                          dtype,
                                          self.sampler.getCurrentIterationNumber() - 1)

                    # update results
                    self.stats.updateResults(dtype, t, self.learner)

            if self.verbose:
                print


    def recomputeStats(self):
        if len(self.dataContainer) == 0:
            self.updateDataContainer(True)

        self.knowledge.clearAlphas()

        for dtype, values in self.dataContainer.items():
            knowledge = {}
            # do the learning
            for t, dataContainer in values.items():
                if dataContainer is not None:
                    for iteration in self.knowledge.getAvailableIterations():
                        grid = self.knowledge.getGrid(self._qoi, iteration)
                        self.learner.grid = grid
                        trainSubset = dataContainer.getTrainDataset()
                        testSubset = None
                        if self.learnWithTest:
                            testSubset = dataContainer.getTestDataset()

                        # -----------------------------------------------------
                        # do the learning
                        if self.verbose:
                            print "learning: t = %g," % t,
                        # learn alpha with corresponding grid
                        self.learner.grid = grid

                        if self.learnWithTest:
                            if self.verbose:
                                print "with test (i=%i, gs=%i)" % (iteration, grid.getSize())
                            self.learner.learnDataWithTest(dataContainer, dtype=dtype)
                        else:
                            if self.verbose:
                                print "(i=%i, gs=%i)" % (iteration, grid.getSize())
                            self.learner.learnDataWithoutTest(dataContainer, dtype=dtype)

                        # -----------------------------------------------------
                        # update the knowledge
                        self.knowledge.update(grid, self.learner.alpha,
                                              self._qoi,
                                              t,
                                              dtype,
                                              iteration)

                        # update stats -> copy
                        self.stats.updateResults(dtype, t, self.learner)

            if self.verbose:
                print
        

    def getParameters(self):
        return self.__params

    def setParameters(self, value):
        self.__params = value

    def getRefinement(self):
        return self._refinement

    def getTestSet(self):
        return self.testSet

    def setTestSet(self, value):
        self.testSet = value

    def setLearnWithTest(self, value):
        self.learnWithTest = value

    def getKnowledgeTypes(self):
        return self.__knowledgeTypes

    def setKnowledgeTypes(self, value):
        self.__knowledgeTypes = value

    def getQoI(self):
        return self._qoi

    def getTimeStepsOfInterest(self):
        return self.__timeStepsOfInterest

    def getRefinement(self):
        return self.refinementManager

    def setRefinement(self, value):
        self.refinementManager = value

    def setQoI(self, value):
        self._qoi = value

    def setTimeStepsOfInterest(self, value):
        self.__timeStepsOfInterest = value

    def getDim(self):
        return self.__params.getDim()

    def getGrid(self):
        return self.sampler.getGrid()

    def setSampler(self, sampler):
        self.sampler = sampler

    def getKnowledge(self):
        return self.knowledge

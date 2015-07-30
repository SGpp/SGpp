from bin.data.DataContainer import DataContainer
from bin.learner.Learner import Learner
from bin.uq.analysis import KnowledgeTypes
from bin.uq.analysis.asgc.ASGCKnowledge import ASGCKnowledge
from pysgpp import DataVector, DataMatrix

from SimulationLearnerSpecification import SimulationLearnerSpecification
import numpy as np
import copy
from bin.uq.operations.sparse_grid import copyGrid, dehierarchize


class SimulationLearner(Learner):

    def __init__(self):
        super(SimulationLearner, self).__init__()
        # specification
        self.specification = SimulationLearnerSpecification()
        self.knowledge = ASGCKnowledge()

        # statistics per knowledge type
        self.trainingOverall = {}
        self.trainAccuracy = {}
        self.trainCount = {}
        self.testAccuracy = {}
        self.testingOverall = {}
        self.testCount = {}

        self.numberPoints = {}
        self.level = {}

        # initialize data container
        self.dataContainer = {}

        # there is one learner per time step
        self._learners = {}

        self._verbose = True

    # ----------------------------------------------------------------
    def __getattr__(self, attr):
        """
        Overrides built-in method if method called is not a object
        method of this Descriptor, most probably it's a method of
        the learner so it tries to call the method
        from our specification
        @param attr: string method name
        @return: method call in specification
        """
        return getattr(self.specification, attr)

    def addLearner(self, learner, t):
        self._learners[t] = copy.copy(learner)

    def getGrid(self):
        return self.grid

    def getKnowledge(self):
        return self.knowledge

    def getSpecification(self):
        return self.specification

    def setSpecification(self, specification):
        self.specification = specification

    def getLearner(self, t=0):
        return self._learners[t]

    # ----------------------------------------------------------------
    def __prepareDataContainer(self, data, name):
        """
        Prepare data for learning
        @param data: dictionary loaded from UQSetting
        @return dictionary {dtype: {t: <DataContainer>}}
        """
        ans = {}
        U = self.getParameters()\
                .activeParams()\
                .getIndependentJointDistribution()
        for dtype in self.getKnowledgeTypes():
            ans[dtype] = {}
            dim = self.grid.getStorage().dim()

            # prepare data container depending on the given knowledge type
            tmp = KnowledgeTypes.transformData(data, U, dtype)

            # load data for all time steps
            for t, values in tmp.items():
                size = len(values)
                mydata = DataMatrix(size, dim)
                sol = DataVector(size)
                for i, (sample, res) in enumerate(values.items()):
                    p = DataVector(sample.getActiveUnit())
                    mydata.setRow(i, p)
                    sol[i] = float(res)
                ans[dtype][t] = DataContainer(mydata, sol, name)
        return ans

    # @profile
    def setDataContainer(self, trainUQSetting, testUQSetting=None):
        """
        Sets the training data container given a UQSetting
        @param trainUQSetting: UQSetting
        @param testUQSetting: UQSetting
        @param dtype: KnowledgeType
        @todo: improve the UQSetting such that it returns just the
        last computed chunk of samples
        """
        # load time steps and quantity of interest
        toi = self.getTimeStepsOfInterest()
        qoi = self.getQoI()
        # load all data points from UQSetting => empty ps
        data = trainUQSetting.getTimeDependentResults(toi, qoi)
        data = self.__prepareDataContainer(data, 'train')
        # set the new data container
        for dtype, values in data.items():
            self.dataContainer[dtype] = {}
            for t, newDataContainer in values.items():
                self.dataContainer[dtype][t] = newDataContainer

        # if there is a test setting given, combine the train and the
        # test data container
        if testUQSetting is not None:
            data = testUQSetting.getTimeDependentResults(toi, qoi)
            data = self.__prepareDataContainer(data, 'test')
            for dtype, values in data.items():
                for t, newDataContainer in values.items():
                    self.dataContainer[dtype][t] = \
                        self.dataContainer[dtype][t].combine(newDataContainer)

    # ----------------------------------------------------------------
    def learnData(self, *args, **kws):
        # learn data
        for dtype, values in self.dataContainer.items():
            knowledge = {}
            print KnowledgeTypes.toString(dtype)
            # do the learning
            for t, dataContainer in values.items():
                print "t = %g, " % t,
                if dataContainer is not None:
                    learner = self._learners[t]
                    # learn data, if there is any available
                    learner.grid = self.getGrid()
                    learner.dataContainer = dataContainer
                    alpha = learner.learnData()

                    # prepare the answer
                    knowledge[t] = copyGrid(learner.grid), DataVector(alpha)
            print
            # update results
            if len(knowledge) > 0:
                self.updateResults(knowledge, dtype)

    def learnDataWithTest(self, dataset=None, *args, **kws):
        # learn data
        for dtype, values in self.dataContainer.items():
            knowledge = {}
            # do the learning
            for t, dataContainer in values.items():
                print "t = %g, " % t,
                if dataContainer is not None:
                    learner = self._learners[t]
                    # learn data, if there is any available
                    learner.grid = self.getGrid()
                    learner.dataContainer = dataContainer
                    alpha = learner.learnDataWithTest(dtype=dtype)

                    # prepare the answer
                    knowledge[t] = copyGrid(learner.grid), DataVector(alpha)
            print
            # update results
            if len(knowledge) > 0:
                self.updateResults(knowledge, dtype)

    def learnDataWithFolding(self, *args, **kws):
        # learn data
        for dtype, values in self.dataContainer.items():
            knowledge = {}

            # do the learning
            for t, dataContainer in values.items():
                if dataContainer is not None:
                    learner = self._learners[t]
                    # learn data, if there is any available
                    learner.grid = self.getGrid()
                    learner.dataContainer = dataContainer
                    alpha = learner.learnDataWithFolding(dtype=dtype)

                    # prepare the answer
                    knowledge[t] = copyGrid(learner.grid), DataVector(alpha)

            # update results
            if len(knowledge) > 0:
                self.updateResults(knowledge, dtype)

    # ----------------------------------------------------------------
    def updateResults(self, knowledge, dtype):
        # run over all learners
        n = len(knowledge)
        for i, (t, (grid, alpha)) in enumerate(knowledge.items()):
            learner = self._learners[t]
            k = -n + i
            # store accuracies of the learning process
            self.trainAccuracy[dtype][t].append(learner.trainAccuracy[k])
            self.trainCount[dtype][t].append(learner.trainCount[k])
            self.trainingOverall[dtype][t].append(np.mean(self.trainAccuracy[dtype][t]))

            # check if there are testing results available
            if len(learner.testAccuracy) == len(learner.trainAccuracy):
                self.testAccuracy[dtype][t].append(learner.testAccuracy[k])
                self.testCount[dtype][t].append(learner.testCount[k])
                self.testingOverall[dtype][t].append(np.mean(self.testAccuracy[dtype][t]))

            # update knowledege
            self.knowledge.update(grid, alpha, self.getQoI(),
                                  t, dtype, self.iteration)
        # update other stats
        gs = grid.getStorage()
        self.level[dtype].append(gs.getMaxLevel())
        self.numberPoints[dtype].append(gs.size())

    def getL2NormError(self):
        """
        calculate L2-norm of error for all learners
        @return: numpy array of L2 norm errors
        """
        return np.array([learner.getL2NormError()
                         for learner in self._learners])

    def getMaxError(self):
        """
        calculate max error for all learners
        @return: numpy array of max error
        """
        return np.array([learner.getMaxError()()
                         for learner in self._learners])

    def getMinError(self):
        """
        calculate min error for all learners
        @return: numpy array of min error
        """
        return np.array([learner.getMinError()()
                         for learner in self._learners])

    # ----------------------------------------------------------------
    def getCollocationNodes(self):
        """
        Create a set of all collocation nodes
        """
        gs = self.grid.getStorage()
        ps = np.ndarray([gs.size(), gs.dim()], dtype='float32')
        p = DataVector(gs.dim())
        for i in xrange(gs.size()):
            gs.get(i).getCoords(p)
            ps[i, :] = p.array()

        return ps

    def refineGrid(self):
        # load the time steps we use for refinement
        # refinets = self.getRefinement().getAdaptTimeWindow()
        refinets = self.getTimeStepsOfInterest()
        oldGridSize = self.getGrid().getSize()
        oldAdmissibleSetSize = self.getRefinement().getAdmissibleSet().getSize()

        # refine
        newCollocationNodes = self.getRefinement().refineGrid(self, refinets)

        # increase counter
        self.iteration += 1

        # print some information
        if self._verbose:
            print "iteration: %i" % self.iteration
            print "old grid size: %i" % oldGridSize
            print "old AS size: %i" % oldAdmissibleSetSize
            print "new collocation nodes: %i" % len(newCollocationNodes)
            print "new grid size:", self.getGrid().getSize()
            print "new AS size: %i" % self.getRefinement()\
                                          .getAdmissibleSet()\
                                          .getSize()

#         fig = plotGrid(self.__grid, self.__knowledge.getAlpha(self.getQoI()),
#                        self.getRefinement().getAdmissibleSetCreator()
#                                            .getAdmissibleSet(),
#                        self.getParameters(), newCollocationNodes)
#         fig.savefig('%i.png' % self._learner.iteration)

        # parse them to a numpy array
        gs = self.grid.getStorage()
        p = DataVector(gs.dim())
        ans = np.ndarray([len(newCollocationNodes), gs.dim()], dtype='float32')
        for i, gp in enumerate(newCollocationNodes):
            gp.getCoords(p)
            ans[i, :] = p.array()

        return ans

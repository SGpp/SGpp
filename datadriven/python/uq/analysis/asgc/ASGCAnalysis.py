#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
#
"""
@file    ASGC.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 23 12:58:31 2013

@brief   Adaptive Sparse Grid Collocation method for UQ

@version  0.1

"""
import os
import json
from anova import HDMRAnalytic
from pysgpp.extensions.datadriven.uq.estimators import MonteCarloStrategy
from pysgpp.extensions.datadriven.uq.operations import (evalSGFunctionMulti,
                               evalSGFunction,
                               isNumerical, discretize)
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF, eval_fullGrid
from pysgpp import (DataVector,
                    DataMatrix)
import numpy as np

from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.tools import writeAlphaARFF, writeGrid
from pysgpp.extensions.datadriven.uq.analysis import Analysis
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation


class ASGCAnalysis(Analysis):
    """
    The ASGC class
    """

    def __init__(self, uqManager, strategy=None):
        Analysis.__init__(self, uqManager.getQoI(),
                          uqManager.getKnowledge().getAvailableTimeSteps(),
                          uqManager.getKnowledge().getAvailableIterations())

        self.__uqManager = uqManager
        self.__params = uqManager.getParameters().activeParams()
        self.__knowledge = uqManager.getKnowledge()
        self.__estimationStrategy = strategy

        # initialize pdf and transformation
        self.__U = self.__params.getIndependentJointDistribution()
        self.__T = self.__params.getJointTransformation()

        # init strategy
        if strategy is None:
            self.__estimationStrategy = MonteCarloStrategy()
        else:
            self.__estimationStrategy = strategy

        # store anova components
        self.__anova = {}

    def getUQManager(self):
        return self.__uqManager

    def getGrid(self):
        qoi = self.__uqManager.getQoI()
        return self.__knowledge.getGrid(qoi)

    def getSurpluses(self):
        qoi = self.__uqManager.getQoI()
        return self.__knowledge.getAlpha(qoi=qoi)

    def setVerbose(self, verbose):
        self.__verbose = verbose

    def eval(self, samples, ts=None, dtype=KnowledgeTypes.SIMPLE):
        ans = {}
        if ts is None:
            ts = self.__uqManager.getTimeStepsOfInterest()

        qoi = self.__uqManager.getQoI()
        grid = self.__knowledge.getGrid(qoi)
        for t in ts:
            alpha = self.__knowledge.getAlpha(qoi, t, dtype)
            ans[t] = evalSGFunction(grid, alpha, samples)

        if len(ts) == 1:
            ans = ans[ts[0]]

        return ans

    # -----------------------------------------------------------------
    # SG probabilistc analysis
    # -----------------------------------------------------------------
    def generateUnitSamples(self, n=10000):
        samples = self.__params.getIndependentJointDistribution().rvs(n)
        trans = self.__params.getJointTransformation()
        return trans.probabilisticToUnitMatrix(samples)

    def estimateDensity(self,
                        ts=[0],
                        n=10000,
                        dtype="kde",
                        samples=None,
                        config={}):
        if samples is None:
            samples = self.generateUnitSamples(n)
        else:
            trans = self.__params.getJointTransformation()
            samples = trans.probabilisticToUnitMatrix(samples)

        time_dependent_values = np.vstack(self.eval(samples, ts=ts))

        if len(ts) == 1:
            return self._estimateDensityByConfig(dtype, time_dependent_values, config)

        ans = {}
        for t, values in time_dependent_values.items():
            ans[t] = self._estimateDensityByConfig(dtype, values, config)

        return ans


    def setEstimationStrategy(self, strategy):
        self.__estimationStrategy.getEstimator().setStrategy(strategy)
        # reset moments
        self.__moments = {}
        self.__moments[self._qoi] = {}

    def computeMean(self, iteration, qoi, t):
        grid, alpha = self.__knowledge.getSparseGridFunction(qoi, t,
                                                             iteration=iteration)
        # do the estimation
        return self.__estimationStrategy.mean(grid, alpha,
                                              self.__U, self.__T)

    def computeVar(self, iteration, qoi, t):
        # compute the mean
        if not self._moments.hasMoment(iteration, qoi, t, 'mean'):
            mean = self.computeMean(iteration, qoi, t)
        else:
            mean = self._moments.getMoment(iteration, qoi, t, 'mean')

        # get the sparse grid function
        grid, alpha = self.__knowledge.getSparseGridFunction(self._qoi, t,
                                                             iteration=iteration)
        # do the estimation
        var = self.__estimationStrategy.var(grid, alpha,
                                            self.__U, self.__T,
                                            mean["value"])

        return {"value": var["value"],
                "err": var["err"] + mean["err"],
                "confidence_interval": var["confidence_interval"]}

    # ----------------------------------------------------------------
    # sensitivity analysis
    # ----------------------------------------------------------------
    def getAnovaDecomposition(self, t=0, iteration=None, *args, **kws):
        # init dictionary
        if iteration is None:
            iteration = self.__knowledge.getIteration()

        if self._qoi not in self.__anova:
            self.__anova[self._qoi] = {}
        if iteration not in self.__anova[self._qoi]:
            self.__anova[self._qoi][iteration] = {}

        # check if the decomposition already exists
        if t in self.__anova[self._qoi][iteration]:
            return self.__anova[self._qoi][iteration][t]

        # determine the anova representation
        grid, alpha = self.__knowledge\
                          .getSparseGridFunction(self._qoi, t=t,
                                                 dtype=KnowledgeTypes.SIMPLE,
                                                 iteration=iteration)
        anova = HDMRAnalytic(grid, alpha, self.__params, *args, **kws)
        anova.setVerbose(self._verbose)
        anova.doDecomposition()

        # store it ...
        self.__anova[self._qoi][iteration][t] = anova

        # ... and return it
        return anova


#     # def histogram(self, t=0, n=1000):
#     #     qoi = self.__specification.getQoI()
#
#     #     grid = self.getGrid()
#     #     alpha = self.__knowledge.getAlpha(qoi, t)
#     #     params = self.__specification.getParameters()
#     #     trans = self.__uqSetting.getSpecification().getTransformation().inv_trans
#
#     #     opEval = createOperationEval(grid)
#     #     for i in xrange(1, n):
#     #         x = trans(params.sample())
#     #         y = opEval.eval(alpha, DataVector(x))
#
#     # ----------------------------------------------------------------
#     # ASGC File Formatter
#     # ----------------------------------------------------------------
#     def setMemento(self, memento):
#         """
#         Restores the state which is saved in the given memento
#
#         Arguments:
#         memento -- the memento object
#
#         Return nothing
#         """
#         self.fromJson(memento)
#
#     def createMemento(self):
#         """
#         Creates a new memento to hold the current state
#
#         Arguments:
#
#         Return a new memento
#         """
#         jsonString = self.toJson()
#         jsonObject = json.JsonReader().read(jsonString)
#         return jsonObject
#
#     @classmethod
#     def fromJson(cls, jsonObject):
#         """
#         Restores the ASGC object from the json object with its
#         attributes.
#
#         Arguments:
#         jsonObject -- json object
#
#         Return the restored ASGC object
#         """
#         setting = ASGC()
#
#         # restore surplusses
#         key = '_ASGC__knowledge'
#         if key in jsonObject:
#             knowledge = ASGCKnowledge.fromJson(jsonObject[key])
#             setting.setKnowledge(knowledge)
#
#         # moments
#         key = '_ASGCSpecification__moments'
#         if key in jsonObject:
#             setting.setMoments(jsonObject[key])
#
#         # restore specification
#         spec = ASGCSpecification()
#
#         # number of moments
#         key = '_ASGCSpecification__k'
#         if key in jsonObject:
#             spec.setNumberOfMoments(int(jsonObject[key]))
#
#         # quantity of interest
#         key = '_ASGCSpecification__qoi'
#         if key in jsonObject:
#             spec.setQoI(jsonObject[key])
#
#         # # adaptive time window
#         # key = '_ASGCSpecification__adaptTimeWindow'
#         # if jsonObject.has_key(key):
#         #     atw = jsonObject[key]
#         #     spec.setAdaptTimeWindow(atw)
#
#         # number of moments
#         key = '_ASGCSpecification__params'
#         if key in jsonObject:
#             params = ParameterSet.fromJson(jsonObject[key])
#             spec.setParameters(params)
#
#         # distribution
#         key = '_ASGCSpecification__distribution'
#         if key in jsonObject:
#             dist = Dist.fromJson(jsonObject[key])
#             spec.setDistribution(dist)
#
#         # distribution
#         key = '_ASGCSpecification__strategy'
#         if key in jsonObject:
#             strategy = EstimationStrategy.fromJson(jsonObject[key])
#             spec.setEstimationStrategy(strategy)
#
#         setting.setSpecification(spec)
#
#         # restore verbose setting
#         key = '_UQSetting__verbose'
#         if key in jsonObject:
#             verbose = jsonObject[key] == 'True'
#             setting.setVerbose(verbose)
#
#         return setting
#
#     def toJson(self):
#         """
#         Returns a string that represents the object
#
#         Arguments:
#
#         Return A string that represents the object
#         """
#         serializationString = '"module" : "' + \
#                               self.__module__ + '",\n'
#         for attrName in dir(self):
#             attrValue = self.__getattribute__(attrName)
#             serializationString += ju.parseAttribute(attrValue, attrName)
#
#         for attrName in dir(self.__specification):
#             attrValue = self.__specification.__getattribute__(attrName)
#             serializationString += ju.parseAttribute(attrValue, attrName)
#
#         s = serializationString.rstrip(",\n")
#
#         # print "j-------------------------------------------"
#         # print "{" + s + "}"
#         # print "j-------------------------------------------"
#
#         return "{" + s + "}"
#
#     def writeErrors(self, radix, ts=None):
#         if not self.__testUQSetting:
#             return
#
#         qoi = self.__specification.getQoI()
#         testData = self.__testUQSetting.getResults(qoi)
#
#         ts = self.getTimeStepsOfInterest()
#
#         # helper variable for pretty printing
#         n = trunc(log(self.__uqSetting.getTimeSetting()['tn'] + 1, 10)) + 2
#
#         for t in ts:
#             if testData and t in testData:
#                 testdataContainer = self.getParameters().getActiveSubset(testData[t])
#
#             alpha = self.__knowledge.getAlpha(qoi, t)
#             grid = self.__grid
#
#             ans = {}
#             for p, val in testdataContainer.items():
#                 ans[p] = evalSGFunction(grid, alpha, DataVector(p)) - val
#
#             # write to file
#             fd = open('%s_t%s_errors.csv' % (radix, str(t).zfill(n)), 'w')
#             for p, val in ans.items():
#                 for pi in p:
#                     fd.write("%f" % pi + ", ")
#                 fd.write("%f" % val + "\n")
#             fd.close()
#
#     def writeGrids(self, filename, ts=None):
#         qoi = self.__specification.getQoI()
#
#         ts = self.getTimeStepsOfInterest()
#         for t in ts:
#             alpha = self.__knowledge.getAlpha(qoi, t)
#             writeCheckpoint("%s.t%g" % (filename, t), self.__grid, alpha)
#
#     def writeStatsRegression(self, filename, ncol=7, t=300):
#         qoi = self.__specification.getQoI()
#         data = DataMatrix(1, ncol)
#         v = DataVector(ncol)
#         v.setAll(0.0)
#
#         v[0] = 0
#         v[1] = self.__grid.getStorage().getMaxLevel()
#         v[2] = self.__grid.getSize()
#         v[3] = 0.0
#
#         # # get alpha and test setting
#         alpha = self.__knowledge.getAlpha(qoi, t)
#         item = self.__testUQSetting.getResults(qoi)[t]
#
#         # # estimate MSE of regression
#         bitmaskActive = [param.isActive() for _, param in self.__specification.getParameters().items()]
#         keys = [[k for i, k in enumerate(key) if bitmaskActive[i]]
#                 for key in item.keys()]
#
#         A = DataMatrix(keys)
#         res = evalSGFunctionMulti(self.__grid, alpha, A)
#         # compute the quadratic error
#         s = sum([(value - fN) ** 2
#                  for value, fN in zip(item.values(), res)])
#         s /= float(len(keys))
#         v[4] = s
#
#         # # estimate expectation value accuracy
#         sg1 = self.E(t)
#         mc1 = sum(item.values()) / float(len(item))
#         v[5] = abs(sg1 - mc1)
#
#         # # estimate variance accuracy
#         sg2 = self.V(t)
#         mc2 = sum([val ** 2 for val in item.values()]) / float(len(item)) - mc1 ** 2
#
#         v[6] = abs(sg2 - mc2)
#
#         data.setRow(0, v)
#
#         writeDataARFF({'filename': filename + ".stats.arff",
#                        'data': data,
#                        'names': ['iteration',
#                                  'level',
#                                  'grid_size',
#                                  'number_of_contributing_points',
#                                  'interpolation_accuracy',
#                                  'expectation_accuracy',
#                                  'variance_accuracy']})
#

# -----------------------------------------------------------------------------
    def computeL2ErrorSurpluses(self, qoi, t, dtype, iteration):
        v1 = np.abs(self.__knowledge.getAlpha(qoi, t, dtype, iteration))
        s = np.sum(v1)

        if iteration > 1:
            v2 = np.abs(self.__knowledge.getAlpha(qoi, t, dtype, iteration - 1))
            s -= np.sum(v2)
        return s

    def computeStats(self, dtype):
        names = ['time',  # 0
                 'iteration',  # 1
                 'level',  # 2
                 'grid_size',  # 3
                 'trainMSE',  # 4
                 'trainL2Error',  # 5
                 'testMSE',  # 6
                 'testL2Error',  # 7
                 'testL1Error',  # 8
                 'testMaxError',  # 9
                 'L2ErrorSurpluses']  # 10

        knowledge = self.__uqManager.getKnowledge()
        ts = knowledge.getAvailableTimeSteps()
        iterations = knowledge.getAvailableIterations()
        nrows = len(ts) * len(iterations)
        ncols = len(names)
        data = DataMatrix(nrows, ncols)
        v = DataVector(ncols)
        v.setAll(0.)
        row = 0

        for t in ts:
            for iteration in iterations:
                v[0] = t
                v[1] = iteration
                v[2] = self.__uqManager.stats.level[dtype][iteration]
                v[3] = self.__uqManager.stats.numberPoints[dtype][iteration]
                v[4] = self.__uqManager.stats.trainMSE[dtype][t][iteration]
                v[5] = self.__uqManager.stats.trainL2Norm[dtype][t][iteration]
                if len(self.__uqManager.stats.testMSE[dtype][t]) == \
                        len(self.__uqManager.stats.trainMSE[dtype][t]):
                    v[6] = self.__uqManager.stats.testMSE[dtype][t][iteration]
                    v[7] = self.__uqManager.stats.testL2Norm[dtype][t][iteration]
                    v[8] = self.__uqManager.stats.testL1Norm[dtype][t][iteration]
                    v[9] = self.__uqManager.stats.testMaxError[dtype][t][iteration]
                v[10] = self.computeL2ErrorSurpluses(self._qoi, t,
                                                     dtype, iteration)
                # write results to matrix
                data.setRow(row, v)
                row += 1
        return {'data': data, 'names': names}

    def writeStats(self, filename):
        for dtype in self.__uqManager.getKnowledgeTypes():
            stats = self.computeStats(dtype)
            suffix = KnowledgeTypes.toString(dtype)
            stats['filename'] = "%s.%s.stats.arff" % (filename, suffix)
            writeDataARFF(stats)

# -----------------------------------------------------------------------------

    def computeMoments(self, iterations=None, ts=None):
        names = ['time',
                 'iteration',
                 'grid_size',
                 'mean',
                 'meanDiscretizationError',
                 'meanConfidenceIntervalBootstrapping_lower',
                 'meanConfidenceIntervalBootstrapping_upper',
                 'var',
                 'varDiscretizationError',
                 'varConfidenceIntervalBootstrapping_lower',
                 'varConfidenceIntervalBootstrapping_upper']
        # parameters
        if ts is None:
            ts = self.__knowledge.getAvailableTimeSteps()
        if iterations is None:
            iterations = self.__knowledge.getAvailableIterations()
        nrows = len(ts) * len(iterations)
        ncols = len(names)
        data = DataMatrix(nrows, ncols)
        v = DataVector(ncols)

        row = 0
        for t in ts:
            for iteration in iterations:
                size = self.__knowledge.getGrid(qoi=self._qoi,
                                                iteration=iteration).getSize()
                mean = self.mean(ts=[t], iterations=[iteration])
                var = self.var(ts=[t], iterations=[iteration])

                v.setAll(0.0)
                v[0] = t
                v[1] = iteration
                v[2] = size
                v[3], v[4] = mean["value"], mean["err"]
                v[5], v[6] = mean["confidence_interval"]
                v[7], v[8] = var["value"], var["err"]
                v[9], v[10] = var["confidence_interval"]

                # write results to matrix
                data.setRow(row, v)
                row += 1

        return {'data': data,
                'names': names}

# -------------------------------------------------------------------------------

    def writeSensitivityValues(self, filename):

        def keymap(key):
            names = self.__uqManager.getParameters().activeParams().getNames()
            ans = [names[i] for i in key]
            return ",".join(ans)

        # parameters
        ts = self.__knowledge.getAvailableTimeSteps()
        gs = self.__knowledge.getGrid(self._qoi).getStorage()

        n = len(ts)
        n1 = gs.getDimension()
        n2 = 2 ** n1 - 1
        data = DataMatrix(n, n1 + n2 + 1)
        names = ['time'] + [None] * (n1 + n2)

        for k, t in enumerate(ts):
            # estimated anova decomposition
            anova = self.getAnovaDecomposition(t=t)
            me = anova.getSobolIndices()

            if len(me) != n2:
                import ipdb; ipdb.set_trace()
            n2 = len(me)
            te = anova.getTotalEffects()
            n1 = len(te)

            v = DataVector(n1 + n2 + 1)
            v.setAll(0.0)
            v[0] = t

            for i, key in enumerate(anova.getSortedPermutations(te.keys())):
                v[i + 1] = te[key]
                if k == 0:
                    names[i + 1] = '"$T_{' + keymap(key) + '}$"'

            for i, key in enumerate(anova.getSortedPermutations(me.keys())):
                v[n1 + i + 1] = me[key]

                if k == 0:
                    names[n1 + 1 + i] = '"$S_{' + keymap(key) + '}$"'

            data.setRow(k, v)

        writeDataARFF({'filename': filename + ".sa.stats.arff",
                       'data': data,
                       'names': names})

#     def sampleGridRandomly(self, filename, n=10000):
#         qoi = self.__specification.getQoI()
#         ts = self.getTimeStepsOfInterest()
#         names = self.__specification.getParameters().getNames()
#         names.append('f_\\mathcal{i}(x)')
#
#         for t in ts:
#             surplus = self.__knowledge.getAlpha(qoi, t)
#
#             # init
#             gs = self.__grid.getStorage()
#             dim = gs.getDimension()
#
#             # do random sampling
#             data = DataMatrix(uniform(0, 1).rvs([n, dim]))
#             res = evalSGFunctionMulti(self.__grid, surplus, data)
#
#             data.transpose()
#             data.appendRow()
#             data.setRow(data.getNrows() - 1, res)
#             data.transpose()
#
#             # -----------------------------------------
#             # write sampled function
#             # -----------------------------------------
#             writeDataARFF({'filename': '%s.t%g.rand_samples.arff' % (filename, t),
#                            'data': data,
#                            'names': names})

    def sampleGrids(self, filename):
        ts = self.__uqManager.getTimeStepsOfInterest()

        names = self.__params.getNames()
        names.append('f_\\mathcal{I}(x)')

        for t in ts:
            grid, surplus = self.__knowledge.getSparseGridFunction(self._qoi, t)

            # init
            gs = grid.getStorage()
            dim = gs.getDimension()

            # -----------------------------------------
            # do full grid sampling of sparse grid function
            # -----------------------------------------
            data = eval_fullGrid(4, dim)
            res = evalSGFunctionMulti(grid, surplus, data)

            data = np.vstack((data.T, res)).T

            # write results
            data_vec = DataMatrix(data)
            writeDataARFF({'filename': "%s.t%f.samples.arff" % (filename, t),
                           'data': data_vec,
                           'names': names})
            del data_vec
            # -----------------------------------------
            # write sparse grid points to file
            # -----------------------------------------
            data = np.ndarray((gs.getSize(), dim))
            x = DataVector(dim)
            for i in xrange(gs.getSize()):
                gp = gs.getPoint(i)
                gs.getCoordinates(gp, x)
                data[i, :] = x.array()

            # write results
            data_vec = DataMatrix(data)
            writeDataARFF({'filename': "%s.t%f.gridpoints.arff" % (filename, t),
                           'data': data_vec,
                           'names': names})
            del data_vec
            # -----------------------------------------
            # write alpha
            # -----------------------------------------
            writeAlphaARFF("%s.t%f.alpha.arff" % (filename, t),
                           surplus)

    def writeCheckpoints(self, filename):
        ts = self.__uqManager.getTimeStepsOfInterest()

        names = self.__params.getNames()
        names.append('f_\\mathcal{I}(x)')

        for iteration in xrange(self.__knowledge.getIteration() + 1):
#             myjson = {"Grid": {"dimNames": ["E", "K_1c", "rho", "n"],
#                                "matrixEntries": ["E", "K_1c", "rho", "n"]},
#                       "Set": {"path": "",
#                               "grids": [],
#                               "alphas": [],
#                               "paramValues": list(ts),
#                               "paramName": "Time"}}
            myjson = {"Grid": {"dimNames": ["phi", "e", "K_L"],
                               "matrixEntries": ["phi", "e", "K_L"]},
                      "Set": {"path": "",
                              "grids": [],
                              "alphas": [],
                              "paramValues": list(ts),
                              "paramName": "Time"}}
            for t in ts:
                grid, surplus = self.__knowledge.getSparseGridFunction(self._qoi, t,
                                                                       iteration=iteration)
                out = "%s.t%f.i%i" % (filename, t, iteration)
                out_grid = "%s.grid" % out
                out_alpha = "%s.alpha.arff" % out
                writeGrid(out_grid, grid)
                writeAlphaARFF(out_alpha, surplus)

                # collect information for json
                myjson["Set"]["grids"].append(os.path.abspath(out_grid))
                myjson["Set"]["alphas"].append(os.path.abspath(out_alpha))

            # write json to file
            fd = open("%s.%i.json" % (filename, iteration), "w")
            json.dump(myjson, fd, indent=2)
            fd.close()

    def computeSurplusesLevelWise(self, t=0, dtype=KnowledgeTypes.SIMPLE,
                                  iteration=None):
        gs = self.__knowledge.getGrid(self._qoi, iteration=iteration).getStorage()
        alpha = self.__knowledge.getAlpha(qoi=self._qoi, t=t, dtype=dtype,
                                          iteration=iteration)

        res = {}
        for i in xrange(gs.getSize()):
            s = gs.getPoint(i).getLevelSum()
            if s not in res:
                res[s] = [alpha[i]]
            else:
                res[s].append(alpha[i])

        return res

    def writeSurplusesLevelWise(self, filename):
        # locate all knowledge types available
        dtypes = self.__uqManager.getKnowledgeTypes()
        names = ['level']
        for dtype in dtypes:
            names.append("surplus_%s" % KnowledgeTypes.toString(dtype))

        ts = self.__knowledge.getAvailableTimeSteps()
        for t in ts:
            # collect all the surpluses classifying them by level sum
            data = {}
            n = 0
            for dtype in dtypes:
                data[dtype] = self.computeSurplusesLevelWise(t, dtype)
                n = sum([len(values) for values in data[dtype].values()])

            A = DataMatrix(n, len(names))
            # add them to a matrix structure
            for i, dtype in enumerate(dtypes):
                k = 0
                for level, surpluses in data[dtype].items():
                    for j, surplus in enumerate(surpluses):
                        A.set(k + j, i + 1, surplus)
                        A.set(k + j, 0, level)
                    k += len(surpluses)

            writeDataARFF({'filename': "%s.t%s.surpluses.arff" % (filename, t),
                           'data': A,
                           'names': names})

#     def writeRefinementEvaluation(self, filename):
#         gs = self.__grid.getStorage()
#         dim = gs.getDimension()
#         stats = self.__uqManager.newPoints
#         numberPoints = self.__uqManager.numberPoints
#         # init
#         names = ['time', 'iteration', 'grid_size'] + \
#                 [x for d in xrange(dim) for x in ['l_%i' % d, 'i_%i' % d]] + \
#                 ['x_%i' % d for d in xrange(dim)] + \
#                 ['surplus']
#
#         data = DataMatrix(0, len(names))
#
#         # run over all available stats
#         for t, values in stats.items():
#             for i, (iteration, ps) in enumerate(values.items()):
#                 for li, (_, coords, surplus) in ps.items():
#                     p = DataVector([t, iteration, numberPoints[i]] +
#                                    [a for l in li for a in l] +
#                                    coords +
#                                    [surplus])
#                     data.appendRow()
#                     data.setRow(data.getNrows() - 1, p)
#
#         # write results
#         writeDataARFF({'filename': "%s.infos.arff" % filename,
#                        'data': data,
#                        'names': names})
#
#     def __str__(self):
#         return self.toJson()

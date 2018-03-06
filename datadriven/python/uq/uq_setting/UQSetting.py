#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    UQSetting.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 13:16:04 2013

@brief It contains the UQSetting class which represents a
non-intrusive UQ setting defined by a preprocession, simulation and
postprocessing step. The preprocessing transforms the initial
parameters and forwards them to the simulation, from which the
results are given to the postprocessing step, where the desired
quantities of interest are computed out of the simulation results.

@version  0.1

"""
from pysgpp.extensions.datadriven.uq.sampler import Sample, SampleType
from math import ceil
from pysgpp import DataVector, DataMatrix
import subprocess
import json
import os
import sys

from pysgpp.extensions.datadriven.tools import readDataARFF
from UQSettingFormatter import UQSettingFormatter
from UQSpecification import UQSpecification
import pysgpp.extensions.datadriven.uq.jsonLib as ju
import numpy as np
import warnings

# for parallelisation
import remote_worker as remote
from multiprocessing import cpu_count


class UQSampleType:
    RAW = 1,
    PREPROCESSED = 2


class UQSetting(object):
    """
    UQSetting class
    """

    def __init__(self):
        """
        Constructor
        """
        self.__specification = UQSpecification()

        self.__stats_samples = {}
        self.__stats_preprocessor = {}
        self.__stats_preprocessor_reverse = {}
        self.__stats_simulation = {}
        self.__stats_postprocessor = {}

        self._verbose = True

        # parallel stuff (taken from mc_berechnung)
        self.parallelprocesses = sum([props['cores'] for props in remote.hosts.values()])
        self.expectedsamplecount = 0  # set to expected number of samples, for parallellisation
        self.__filesuffix = 0  # incremented for each process
        self.children = {}  # running worker processes
        self.files = []  # result files written
        self.lastid = 0

        # speed up result loading
        self.__dictResults = {}

    def __getattr__(self, attr):
        """
        Overrides built-in method if method called is not a object
        method of this Descriptor, most probably it's a method of the
        specification so it tries to call the method from UQSpecification.
        @param attr: String method name
        """
        return getattr(self.__specification, attr)

    def findEquivalent(self, sample, stats):
        # get sample in unit space
        p = tuple(sample.getExpandedUnit())
        found = p in stats
        if not found:
            min_diff = 10.
            # get parameter in probabilistic space
            x = tuple(sample.getExpandedProbabilistic())
            # DAMN IT => this is shit!!!!
            keys = stats.keys()
            if self._verbose:
                print "search for equivalent for %s" % (p,)
            j = 0
            while not found and j < len(keys):
                g = keys[j]
                j += 1
                diff1 = [abs(gi - pi) / pi for gi, pi in zip(g, p) if abs(pi) > 0]
                diff2 = [abs(gi - pi) / gi for gi, pi in zip(g, p) if abs(gi) > 0]
                diff = sum(diff1) + sum(diff2)
                min_diff = min(min_diff, diff)
                if diff < 1e-6:
                    found = True
                    old_q = self.__stats_preprocessor[g]
                    if self._verbose:
                        print g, old_q

                    # save old items
                    q = self.getPreprocessor().unitToProbabilistic(x)

                    diff1 = [abs(gi - pi) / pi for gi, pi in zip(q, old_q) if abs(pi) > 0]
                    diff2 = [abs(gi - pi) / gi for gi, pi in zip(q, old_q) if abs(gi) > 0]
                    diff = sum(diff1) + sum(diff2)

                    simulation = self.__stats_simulation[old_q].copy()
                    post = self.__stats_postprocessor[old_q].copy()

                    # delete old ones if
                    if g in self.__stats_samples:
                        del self.__stats_samples[g]
                    del self.__stats_preprocessor[g]
                    del self.__stats_preprocessor_reverse[old_q]
                    del self.__stats_simulation[old_q]
                    del self.__stats_postprocessor[old_q]

                    # insert old items with new key
                    self.__stats_samples[p] = sample
                    self.__stats_preprocessor[p] = q
                    self.__stats_preprocessor_reverse[q] = p
                    self.__stats_simulation[q] = simulation
                    self.__stats_postprocessor[q] = post

            if self._verbose:
                if found:
                    print "found equivalent",
                    print min_diff
                    print len(self.__stats_samples), \
                        len(self.__stats_preprocessor), \
                        len(self.__stats_preprocessor_reverse), \
                        len(self.__stats_simulation), \
                        len(self.__stats_postprocessor)
                else:
                    print "no equivalent found"
        return found

    def __preprocessing(self, sample, *args, **kws):
        """
        Transforms a parameter taken from the unit hyper cube into the
        computational range of the problem at hand.
        @param sample: tuple input parameter in [0, 1]^d
        @return: return input parameter for the simulation run
        """
        # get sample in unit space
        p = tuple(sample.getExpandedUnit())
        # get parameter in probabilistic space
        x = tuple(sample.getExpandedProbabilistic())
        preprocessor = self.getPreprocessor()

        # transform the data point
        if preprocessor:
            if p not in self.__stats_preprocessor and \
                    not self.findEquivalent(sample, self.__stats_preprocessor):
                try:
                    t0 = self.getStartTime()
                    tn = self.getEndTime()
                    dt = self.getTimeStep()

                    # do the pre-processing
                    q = preprocessor.unitToProbabilistic(x, *args,
                                                         t0=t0, tn=tn, dt=dt,
                                                         **kws)

                    if self._verbose:
                        print "Apply pre-processing:",
                except:
                    raise
            else:
                q = self.__stats_preprocessor[p]

                if self._verbose:
                    print "Restore pre-processing:",
        else:
            q = x

            if self._verbose:
                print "No pre-processor defined:",

        if self._verbose:
            print "%s -> %s -> %s -> %s" % (tuple(sample.getActiveUnit()), p, x, q)

        # check if pre-processed data has the right dimension
        if len(self.__stats_preprocessor) > 0:
            val = self.__stats_preprocessor.itervalues().next()
            if type(val) != type(q):
                raise TypeError('The pre-processor has changed since',
                                'the last run')
            try:
                a, b = len(val), len(q)
            except TypeError:
                a, b = 0, 0

            if a != b:
                raise TypeError('The pre-processor has changed since',
                                'the last run')

        if p not in self.__stats_preprocessor:
            # store results
            self.__stats_samples[p] = sample
            self.__stats_preprocessor[p] = q

            if q in self.__stats_preprocessor_reverse:
                # the pre processor is surjective => something has
                # changed in the parameter configuration. We delete
                # the old entry in the preprocessor stats
                w = self.__stats_preprocessor_reverse[q]
                del self.__stats_preprocessor[w]
                # raise AttributeError('The transformation function is surjective')

            self.__stats_preprocessor_reverse[q] = p

        assert len(self.__stats_preprocessor) == \
            len(self.__stats_preprocessor_reverse)

        return q

    def __eval(self, q, *args, **kws):
        """
        Run the simulation with the given parameter q
        @param q: tuple input parameter for the simulation run
        @return: simulation result
        """
        simulation = self.getSimulation()
        if simulation:
            # check whether the simulation run has already been
            # performed with the current parameter set
            if q not in self.__stats_simulation:
                if self._verbose:
                    print "Run simulation..."

                try:
                    t0 = self.getStartTime()
                    tn = self.getEndTime()
                    dt = self.getTimeStep()

                    A = simulation(q, *args,
                                   t0=t0, tn=tn, dt=dt,
                                   **kws)

                    self.__stats_simulation[q] = A
                except:
                    raise
            else:
                A = self.__stats_simulation[q]
                if self._verbose:
                    print "Restore simulation result..."
        else:
            raise Exception('Simulation is missing')

        if not len(self.__stats_preprocessor) == len(self.__stats_simulation):
            print q, len(self.__stats_preprocessor), \
                len(self.__stats_simulation)
            # import ipdb; ipdb.set_trace()
            # assert len(self.__stats_preprocessor) == len(self.__stats_simulation)

        return A

    def __postprocessing(self, A, q, *args, **kws):
        """
        Calculates quantities of interest depending on the output of
        the simulation.
        @param A: result of the simulation
        @param q: input parameter used for keying
        @return: quantities of interest defined by the postprocessing
        function
        """
        postprocessor = self.getPostprocessor()
        if postprocessor:
            # check whether the postprocessing part has already been
            # performed with the current parameter set
            if q not in self.__stats_postprocessor:
                if self._verbose:
                    print "Apply post processing..."
                try:
                    t0 = self.getStartTime()
                    tn = self.getEndTime()
                    dt = self.getTimeStep()
                    B = postprocessor(A, *args,
                                      q=q,
                                      t0=t0, tn=tn, dt=dt,
                                      **kws)
                except:
                    raise
            else:
                B = self.__stats_postprocessor[q]

                if self._verbose:
                    print "Restore postprocessing data"
        else:
            B = A
            if self._verbose:
                print "No postprocessor defined"

        # check if post-processed data has the right dimension
        if len(self.__stats_postprocessor) > 0:
            val = self.__stats_postprocessor.itervalues().next()
            if type(val) != type(B):
                raise TypeError('The post-processor has changed since the last run')

        if q not in self.__stats_postprocessor:
            # store result
            self.__stats_postprocessor[q] = B

        return B

    def writeToFile(self, filename=None):
        # save setting to file
        if not filename:
            filename = self.getFilename()
        if not filename:
            if self._verbose:
                print "Filename not set, memento not written."
            return
        if self._verbose:
            print "Write memento to file..."
        m = self.createMemento()
        UQSettingFormatter().serializeToFile(m, filename)

    def runSimulations(self):
        """
        Run the simulations for all available pre-processor stats
        """
        acc = []
        self.__stats_simulation = {}
        for i, (_, q) in enumerate(self.__stats_preprocessor.items()):
            print "run simulation %i/%i" % (i + 1, len(self.__stats_preprocessor.items()))
            ans = self.__eval(q)
        return acc

    def runPostprocessor(self):
        """
        Run the post-processing for all available simulation results.
        """
        acc = []
        self.__stats_postprocessor = {}
        for i, (q, A) in enumerate(self.__stats_simulation.items()):
            print "run postprocessing %i/%i" % (i + 1, len(self.__stats_simulation.items()))
            self.__postprocessing(A, q)
            # if len(ans['time']) == 0:
            #     # invalid setting remove it
            #     p = self.__stats_preprocessor_reverse[q]
            #     del self.__stats_preprocessor_reverse[q]
            #     del self.__stats_preprocessor[p]
            #     del self.__stats_simulation[q]
            #     del self.__stats_postprocessor[q]
            #     acc.append((q, A['massflux']))
        return acc

    def runSamples(self, samples, dist=False, *args, **kws):
        # remove those samples from sample list which have already
        # been evaluated
        samples.removeSet(self.__stats_samples.values())
        if len(samples) == 0:
            return 0

        # do the distributed run
        n = self.getSize("simulation")
        if dist is True:
            self.runSamples_dist(samples)
            self.waitForResults()
            self.loadResults()
        else:
            self.runSamples_withoutDistribution(samples)

        return self.getSize("simulation") - n

    def runSamples_dist(self, samples, starti=0, *args, **kws):
        """
        Performs run() for the given samples.
        Parallelizes automatically, so use this if possible.
        Sampling will take place in the background, use waitForResults()
        at the end.
        """
        if len(samples) == 0:
            return

        # split into smaller chunks
        jobsize = float(self.expectedsamplecount) / float(self.parallelprocesses)
        if jobsize < 1:
            jobsize = 1.0
        if len(samples) > jobsize:
            njobs = ceil(len(samples) / jobsize)  # round up
            nsamples = int(ceil(len(samples) / njobs))
            njobs = int(njobs)
            print "jobconfig:", njobs, jobsize, nsamples, len(samples)

            # split into enough chunks to be below jobsize
            for i in range(0, njobs - 1):
                if (nsamples * (i + 1)) >= len(samples):
                    print i, "are enough"
                    break
                print "job", nsamples * i, len(samples[nsamples * i:nsamples * (i + 1)])
                self.runSamples_dist(samples[nsamples * i:nsamples * (i + 1)],
                                     starti=nsamples * (njobs - 1), *args, **kws)
            self.runSamples_dist(samples[nsamples * (njobs - 1):],
                                 starti=nsamples * (njobs - 1), *args, **kws)
            return

        if 'starti' in kws:
            starti = kws['starti']
        else:
            starti = self.__nextId(len(samples))

        # choose host and do load balancing
        host = ""
        while host == "":
            # look for returned children
            # print "self.children: ", self.children
            children_to_remove = []  # collect the finished children here
            for child_pid, props in self.children.iteritems():
                (pid, res) = os.waitpid(child_pid, os.WNOHANG)
                if pid != 0:
                    # the child already returned
                    if res == 0:
                        self.files.append(self.children[pid][1])

                        if self.children[pid][2] != '':
                            remote.free_host(self.children[pid][2])

                        # del self.children[pid]
                        children_to_remove.append(pid)
                        print "child finished, pid: %d" % (pid)
                    else:
                        starti, f, host = self.children[pid]

                        if host != '':
                            remote.free_host(host)

                        print "WARNING: child %d crashed, File %s, starting from index %d on %s" % (pid, f, starti, host)
                        # del self.children[pid]
                        children_to_remove.append(pid)
            # remove the collected children
            for pid in children_to_remove:
                del self.children[pid]

            host = remote.choose_host()  # try to choose host
            # this will return "" if there are no free cores

        # choose unique working directory on host
        cmd = ["ssh", host,
               "python2 -c 'from pysgpp.extensions.datadriven.uq.uq_setting.choose_working_dir import choose_working_directory; choose_working_directory(\"%s\")'"
               % remote.hosts[host]['scratch'] ]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p.wait()
        scratch_path = p.stdout.next().strip()

        filename = self.getFilename().replace(".gz", "." + str(self.__filesuffix) + ".gz")
        res = os.fork()

        if res == 0:
            # in the child, set a new file and do work
            remote.do_dist_run(host, self, filename, samples, starti, scratch_path, self.__filesuffix)
        elif res > 0:
            # in the parent, clean up and return
            self.children[res] = (starti, filename, host)
            # print "incremented suffix", map(lambda it: (it[0], it[1][0], it[1][2], self.children.iteritems())
            self.__filesuffix = self.__filesuffix + 1
            return
        else:
            print "Error while forking:", starti

#    def computeResultsOfChild(self, pid):
#        """
#        If a child returned, do organization.
#        """


    def waitForResults(self):
        """
        Wait for all forked workers to finish working,
        then load all result files into one SamplingResult.
        """
        # this does not actually belong here, but here it will run in parallel
        self.writeToFile()

        while len(self.children) > 0:
            print "Now waiting for child to return..."
            pid, res = os.wait()
            if res == 0:
                self.files.append(self.children[pid][1])

                if self.children[pid][2] != '':
                    remote.free_host(self.children[pid][2])

                del self.children[pid]
                print "child finished, %d remaining" % len(self.children)
            else:
                starti, f, host = self.children[pid]

                if host != '':
                    remote.free_host(host)

                print "WARNING: child %d crashed, File %s, starting from index %d on %s" % (pid, f, starti, host)
                del self.children[pid]

    def loadResults(self):
        """
        Load Results from files, if necessary
        and merge them with this uq setting
        """
        if self.files != []:
            if len(self.children) > 0:
                self.waitForResults()
            for f in self.files:
                self.loadfile(f)
            self.files = []

    def loadfile(self, filename):
        wd = os.getcwd()
        pwd = wd.split('/')
        pwd[2] = remote.username
        pwd = "/".join(pwd)
        print "scp helium:%s %s" % (os.path.join(pwd, filename), wd)
        subprocess.call("scp helium:%s %s" % (os.path.join(pwd, filename), wd),
                        shell=True)

        m = UQSettingFormatter().deserializeFromFile(filename)
        uq = UQSetting.fromJson(m)
        self.mergeStats(uq)

    def runSamples_withoutDistribution(self, samples, *args, **kws):
        """
        Run for multiple samples
        @param samples: list of samples
        """
        cnt = 0
        for i, sample in enumerate(samples):
            if self._verbose:
                print "-" * 60
                print "Run: %i/%i (%i)" % (i + 1, len(samples), self.getSize())

            ret = self.run(sample, writeToFile=self.getSaveAfterEachRun(cnt),
                           *args, **kws)

            # check if run was successful
            if ret == 1:
                print "Warning: invalid sample %s" % sample
            else:
                cnt += 1

#             if self._verbose:
#                 sys.out.flush()

        if self._verbose:
            print "-" * 60

    def run(self, sample, writeToFile=False, *args, **kws):
        """
        Trigger the pipeline from
         1) Transformation p \in [0, 1]^d -> \Gamma*
         2) Pre-processing
         3) Simulation run
         4) Post-processing
        @param p: tuple input parameter in [0, 1]^d
        @return: simulation results
        """
        try:
            n1 = len(self.getSimulationStats())
    
            # pre processing
            q = self.__preprocessing(sample, *args, **kws)
    
            # do simulation run
            A = self.__eval(q, *args, **kws)
    
            n2 = len(self.getSimulationStats())
    
            # serialize the UQSetting to file. As this is the expensive
            # part, this step is important if there is any error in the post
            # processing part
            if n2 > n1 and writeToFile:
                self .writeToFile()
    
            n1 = len(self.getPostprocessorStats())
    
            self.__postprocessing(A, q, *args, **kws)
    
            n2 = len(self.getPostprocessorStats())
            # serialize once more after the post processing has been done
            if n2 > n1 and writeToFile:
                self.writeToFile()
    
            return 0
        except Exception as exception:
            print str(exception)
            return 1

    def sanityCheck(self):
        """
        Check for consistency of the UQ setting.
        """
        ans = 0
        for p, val in self.__stats_preprocessor.items():
            if val not in self.__stats_preprocessor_reverse:
                ans += 1
                print "pre %s" % (p,)
            if val not in self.__stats_simulation:
                ans += 1
                print "sim %s" % (p,)
            if val not in self.__stats_postprocessor:
                ans += 1
                print "post %s" % (p,)
        return ans

    def split(self, n=500):
        """
        Split the UQ setting in two disjunct parts
        @param n: number of elements of the new UQ setting
        """
        ans = UQSetting()
        ans.setSpecification(self.__specification)
        ans.setVerbose(self._verbose)
        npre = {}
        npreRev = {}
        nsim = {}
        npost = {}

        pre = self.__stats_preprocessor.items()
        ixs = range(len(self.__stats_preprocessor))[:n]
        for i in ixs:
            # copy to new object
            p, val = pre[i]
            npre[p] = val
            npreRev[val] = p
            nsim[val] = self.__stats_simulation[val]
            npost[val] = self.__stats_postprocessor[val]

            # delete from old one
            del self.__stats_preprocessor[p]
            del self.__stats_preprocessor_reverse[val]
            del self.__stats_simulation[val]
            del self.__stats_postprocessor[val]

        ans.setPreprocessorStats(npre)
        ans.setPreprocessorStatsReverse(npreRev)
        ans.setSimulationStats(nsim)
        ans.setPostprocessorStats(npost)

        # # change filenames
        # ans.getSpecification().setFilename('1-' + ans.getSpecification().getFilename())
        # self.getSpecification().setFilename('2-' + self.getSpecification().getFilename())

        return ans

    def getAvailableQoI(self):
        """
        Get the availbale quantities of interest.
        @return: list of strings containing the available quantities of
        interest defined by the post-processor function
        """
        if self.__stats_postprocessor:
            # get some arbitrary result from which the keys are the
            # available quantities of interest
            values = self.__stats_postprocessor.values()
            qois = []
            i = 0
            while i < len(values) and len(qois) == 0:
                qois = values[i].keys()
                i += 1
            return qois
        else:
            return []

    def getResult(self, sample, ts=[0], qoi='_'):
        """
        get the result for one given sample
        @param sample: Sample
        @param ts: list of numerics
        @param qoi: string
        @result: numpy array with the scalar results per time step
        """
        p = tuple(sample.getExpandedUnit())

        if p not in self.__stats_preprocessor:
            raise AttributeError('there are no results available for %s' % (p,))

        q = self.__stats_preprocessor[p]

        if q not in self.__stats_postprocessor:
            return None

        if qoi not in self.getAvailableQoI():
            raise AttributeError('quantity of interest "%s" does not exist. \
                available are %s' % (qoi, self.getAvailableQoI()))

        # load available time setting for parameter q
        stats_ts = self.getTimeSetting(q)['t']

        # check if there is a result available
        results = self.__stats_postprocessor[q][qoi]

        # assert that both, results and time stamps have the same length
        if len(results) != len(stats_ts) and len(results) != len(ts):
            raise AttributeError('number of simulation runs (%i) differ from \
                                  number of time steps (%i) for %s' %
                                  (len(results), len(stats_ts), q))

        # make sure that there are results
        if len(stats_ts) == 0:
            return None

        # search for explicitly given results
        ans = np.ndarray([len(ts)], dtype='float')

        # just return the value if there is just one available
        if len(stats_ts) == 1:
            ans[0] = results[0]
            return ans

        # sort it to search it faster
        stats_ixs = np.argsort(stats_ts)
        stats_ts = [stats_ts[ix] for ix in stats_ixs]
        if len(results) == len(stats_ts):
            results = [results[ix] for ix in stats_ixs]

        # split ts in inter/extrapolating time steps
        inter_ixs = [i for i, t in enumerate(ts)
                     if stats_ts[0] <= t <= stats_ts[-1]]
        extraUpper_ixs = [i for i, t in enumerate(ts) if t > stats_ts[-1]]
        extraLower_ixs = [i for i, t in enumerate(ts) if t < stats_ts[0]]

        # fill lower extrapolation with zeros
        for ix in extraLower_ixs:
            ans[ix] = 0.

        # fill upper extrapolation if a steady state is defined
        if len(extraUpper_ixs) > 0:
            if self.reachesSteadyState():
                # do constant extrapolation
                for i in extraUpper_ixs:
                    ans[i] = results[-1]
            else:
                # error: extrapolation needed, but not specified
                raise AttributeError('extrapolation needed (%g < %g), but\
                                     no steady state reached' %
                                     (stats_ts[-1], ts[extraUpper_ixs[0]]))

        # load the results for which the results is directly available
        interpolated_ts = {}
        for i in inter_ixs:
            t = ts[i]
            ix = np.searchsorted(stats_ts, t)
            if ix < len(stats_ts) and abs(t - stats_ts[ix]) < 1e-8:
                ans[i] = results[ix]
            else:
                interpolated_ts[i] = t

        # if we reach this area, no explicit result was found => check if
        # there is an interpolation method specified
        if len(interpolated_ts) > 0 and len(results) > 1 and \
                self.hasInterpolationFunction():
            # sort the input for interpolation
            interpolant = self.getInterpolationFunction(p, stats_ts, results)

            # evaluate interpolation function at each missing time step
            for i, t in interpolated_ts.items():
                # load them now by interpolation
                ans[i] = interpolant(t)

        return ans

    def hasResult(self, sample):
        p = tuple(sample.getExpandedUnit())
        if p in self.__stats_preprocessor:
            q = self.__stats_preprocessor[p]
            if q in self.__stats_simulation and \
                    q in self.__stats_postprocessor:
                return True
        return False

    def getTimeDependentResults(self, ts, qoi='_', ps=None,
                                sampleType=UQSampleType.RAW):
        """
        Collects the simulation results for the given time steps and a certain
        quantity of interest. If just a selection of parameters is needed, one
        can specify it using the ps parameter.
        @param ts: list of time steps
        @param ps: Samples, selection of parameters
        @param qoi: quantity of interest
        @param ps: list of samples to be loaded

        @return: dictionary {<time step>: {<Sample>: value}}
        """
        keyResults = (tuple(ts), qoi)
        if ps is None and keyResults in self.__dictResults and \
                np.all([len(self.__dictResults[keyResults][t]) == self.getSize()
                        for t in ts]):
            return self.__dictResults[keyResults]
            
        if qoi not in self.getAvailableQoI():
            raise AttributeError(('the quantity of interest "%s" does not ' +
                                  'exist. There are "%s" available.') %
                                 (qoi, self.getAvailableQoI()))

        if ps is None:
            # collect all available results -> just the ones which have
            # an entry in the post processing table
            ps = [None] * len(self.__stats_postprocessor)
            for i, q in enumerate(self.__stats_postprocessor):
                p = self.__stats_preprocessor_reverse[q]
                if p not in self.__stats_samples:
                    raise AttributeError('Sample is missing for %s' % (p,))

                ps[i] = self.__stats_samples[p]

        if self._verbose:
            print "loading results %i x %i = %i" % \
                (len(ts), len(ps), len(ts) * len(ps))

        B = {}
        # init B
        for t in ts:
            B[t] = {}

        # run over all parameters
        for k, p in enumerate(ps):
            if self._verbose:
                print " %i/%i \r" % (k + 1, len(ps)),

            res = self.getResult(p, ts, qoi)
            
            # select the key
            if sampleType == UQSampleType.PREPROCESSED:
                key = self.__stats_preprocessor[tuple(p.getExpandedUnit())]
            else:
                # sampleType == UQSampleType.RAW
                key = p
            
            if res is not None:
                # run over all time steps
                for i, t in enumerate(ts):
                    B[t][key] = res[i]

        if self._verbose:
            print

        # store result
        self.__dictResults[keyResults] = B

        return B

    def getResults(self, qoi="_", sampleType=UQSampleType.RAW):
        """
        Collects the simulation results assuming that ts = [0].
        If just a selection of parameters is needed, one
        can specify it using the ps parameter.
        @param sampleType:
        @return: dictionary {0: {<Sample>: value}}
        """
        # collect samples
        ans = {}

        # collect all available results -> just the ones which have
        # an entry in the post processing table
        ps = [None] * len(self.__stats_postprocessor)
        for i, (q, value) in enumerate(self.__stats_postprocessor.items()):
            p = self.__stats_preprocessor_reverse[q]
            if p not in self.__stats_samples:
                raise AttributeError('Sample is missing for %s' % (p,))

            # select the key
            sample = self.__stats_samples[p]
            if sampleType == UQSampleType.PREPROCESSED:
                key = self.__stats_preprocessor[tuple(self.__stats_samples[p].getExpandedUnit())]
            else:
                # sampleType == UQSampleType.RAW
                key = p

            ans[key] = value[qoi]
        
        return {0: ans}

    def getSamplesStats(self):
        return self.__stats_samples

    def setSamplesStats(self, value):
        self.__stats_samples = value

    def getPreprocessorStats(self):
        """
        Get the pre-processor results

        @return: dictionary [<unit parameter>] -> [<paramter>]
        containing the pre-processor results
        """
        return self.__stats_preprocessor

    def getPreprocessorStatsReverse(self):
        """
        Get the reversed pre-processor results
        @return: dictionary [<parameter>] -> [<unit-parameter>]
        of reversed pre-processor stats
        """
        return self.__stats_preprocessor_reverse

    def getSimulationStats(self):
        """
        Get the simulation results
        @return: dictionary [<paramter>] -> simulation result
        """
        return self.__stats_simulation

    def getPostprocessorStats(self):
        """
        Get post-processor results
        @return: dictionary [<parameter>] -> [<qoi>][...]
        containing the post-processor results
        """
        return self.__stats_postprocessor

    def setPreprocessorStats(self, stats):
        """
        Set pre-processor stats
        @param stats: dictionary
        """
        self.__stats_preprocessor = stats

    def setPreprocessorStatsReverse(self, stats):
        """
        Set reversed pre-processor stats
        @param stats: dictionary
        """
        self.__stats_preprocessor_reverse = stats

    def setSimulationStats(self, stats):
        """
        Set simulation results
        @param stats: dictionary
        """
        self.__stats_simulation = stats.copy()

    def setPostprocessorStats(self, stats):
        """
        Set post-processor results
        @param stats: dictionary
        """
        self.__stats_postprocessor = stats

    def setVerbose(self, verbose):
        """
        Set verbose level
        @param verbose: bool
        """
        self._verbose = verbose

    def setLastId(self, lastid):
        self.lastid = lastid

    def __nextId(self, number=1):
        myid = self.lastid + number
        self.lastid = myid
        return myid

    def getDim(self):
        """
        @return: stochastic dimension
        """
        if self.__specification:
            return self.getDim()

    def __len__(self):
        """
        @return: number of simulation results
        """
        return len(self.__stats_simulation)

    def setSpecification(self, specification):
        """
        Set the specification
        @param specification: UQSpecification object
        """
        self.__specification = specification

    def getSpecification(self):
        """
        @return: The UQSpecification object
        """
        return self.__specification

    def hasTimeSetting(self):
        """
        @return: bool defining whether there is a specified time setting or not
        """
        ts = self.getTimeSetting()
        return ts['t0'] < ts['tn']

    def getTimeSetting(self, q=None):
        """
        Get time setting consisting of 't0', 'tn', and 'dt'
        @return: dictionary containing the time setting for the simulation
        """
        t0 = self.getStartTime()
        tn = self.getEndTime()
        dt = self.getTimeStep()
        n = 0

        # get time steps by point
        t = {}
        if q is not None and q in self.__stats_postprocessor and \
                'time' in self.__stats_postprocessor[q]:
            t[q] = self.__stats_postprocessor[q]['time']
            n = len(t)
        else:
            # {p : {qoi: [...]}}
            for p, val in self.__stats_postprocessor.items():
                # search for specific time setting
                if 'time' in val:
                    t[p] = val['time']

                # maximal number of time steps available
                n = max(n, len(val.values()[0]))

        if len(t) > 0:
            # flatten all the available time steps
            ts = t.values()
            if len(ts) == 0:
                return {'t': []}
            t0 = float('inf')
            tn = float('-inf')
            for ti in ts:
                if len(ti) == 0:
                    return {'t': []}
                t0 = min(t0, min(ti))
                tn = max(tn, max(ti))

            dt = None
            if q is not None:
                t = t[q]
        elif dt > 0 and t0 >= 0 and t0 < tn:
            # valid time setting found
            pass
        elif dt > 0 and t0 > -1 and tn == -1:
            # start time and time step are given, but end time is
            # missing => extract it from simulation results if there
            # are any
            tn = (n - 1) * dt
        elif dt == -1 and t0 == -1 and tn == -1:
            # no time setting => use standard setting
            t0 = 0
            tn = n - 1
            dt = 1
        else:
            raise Exception('Invalid time setting found \
                            (%s, %s, %s)' % (t0, tn, dt))

        if len(t) == 0:
            # create t setting
            n = int(ceil((tn - t0) / dt + 1))
            t = np.linspace(t0, tn, n, endpoint=True)

        # return for given time steps
        return {'t0': t0, 'tn': tn, 'dt': dt, 'n': n, 't': t}

    def getVerbose(self):
        """
        @return: verbosity level
        """
        return self._verbose

    def getSize(self, item="postprocessor"):
        """
        @oaram item: string identifying the stage
        @return: int number of simulation results
        """
        if item == "preprocessor":
            return len(self.__stats_preprocessor)
        elif item == "simulation":
            return len(self.__stats_simulation)
        elif item == "postprocessor":
            return len(self.__stats_postprocessor)
        elif item == "samples":
            return len(self.__stats_samples)
        else:
            raise AttributeError("item attribute '%s' unknown" % item)

    def mergeStats(self, newSetting):
        """
        Merge statistics in this setting and some other.
        @param newSetting: UQSetting to be merged
        """
        statsSamples = newSetting.getSamplesStats()
        if statsSamples and self.__stats_samples is not None:
            for key, val in statsSamples.items():
                if key not in self.__stats_samples:
                    self.__stats_samples[key] = val
                else:
                    raise AttributeError('overlapping settings => cannot \
                                          merge them automatically')

        statsPreprocessor = newSetting.getPreprocessorStats()
        if statsPreprocessor and self.__stats_preprocessor is not None:
            for key, val in statsPreprocessor.items():
                if key not in self.__stats_preprocessor:
                    self.__stats_preprocessor[key] = val
                else:
                    raise AttributeError('overlapping settings => cannot \
                                          merge them automatically')

        statsPreprocessorReverse = newSetting.getPreprocessorStatsReverse()
        if statsPreprocessorReverse and self.__stats_preprocessor_reverse is not None:
            for key, val in statsPreprocessorReverse.items():
                if key not in self.__stats_preprocessor_reverse:
                    self.__stats_preprocessor_reverse[key] = val
                else:
                    raise AttributeError('overlapping settings => cannot \
                                          merge them automatically')

        statsSimulation = newSetting.getSimulationStats()
        if statsSimulation and self.__stats_simulation is not None:
            for key, val in statsSimulation.items():
                if key not in self.__stats_simulation:
                    self.__stats_simulation[key] = val
                else:
                    raise AttributeError('overlapping settings => cannot \
                                          merge them automatically')

        statsPostprocessor = newSetting.getPostprocessorStats()
        if statsPostprocessor and self.__stats_postprocessor is not None:
            for key, val in statsPostprocessor.items():
                if key not in self.__stats_postprocessor:
                    self.__stats_postprocessor[key] = val
                else:
                    raise AttributeError('overlapping settings => cannot \
                                          merge them automatically')

        self.lastid = max(newSetting.lastid, self.lastid)

    def merge(self, newSetting):
        """
        Merges this UQSetting with the one in newSetting
        @param newSetting: UQSetting to be merged
        """
        # merge UQSetting
        statsSamples = newSetting.getSamplesStats()
        if statsSamples and not self.__stats_samples:
            self.__stats_samples = statsSamples

        statsPreprocessor = newSetting.getPreprocessorStats()
        if statsPreprocessor is not None and self.__stats_preprocessor is None:
            self.__stats_preprocessor = statsPreprocessor

        statsPreprocessorReverse = newSetting.getPreprocessorStatsReverse()
        if statsPreprocessorReverse is not None and self.__stats_preprocessor_reverse is None:
            self.__stats_preprocessor_reverse = statsPreprocessorReverse

        statsSimulation = newSetting.getSimulationStats()
        if statsSimulation is not None and self.__stats_simulation is None:
            self.__stats_simulation = statsSimulation

        statsPostprocessor = newSetting.getPostprocessorStats()
        if statsPostprocessor is not None and self.__stats_postprocessor is None:
            self.__stats_postprocessor = statsPostprocessor

        # merge specification
        specification = newSetting.getSpecification()

        if specification:
            t0 = specification.getStartTime()
            if t0 != -1 and self.getStartTime() == -1:
                self.setStartTime(t0)

            tn = specification.getEndTime()
            if tn != -1 and self.getEndTime() == -1:
                self.setEndTime(tn)

            dt = specification.getTimeStep()
            if dt != -1 and self.getTimeStep() == -1:
                self.setTimeStep(dt)

        # # check if parameters have changed and adjust the statistics
        # # accordingly. Just a change of uncertain to deterministic and
        # # vice versa is supported.

        # oldParams = specification.getParameters()
        # newParams = self.getParameters()

        # if oldParams and newParams and oldParams != newParams:
        #     oldTrans = specification.getTransformation()
        #     newTrans = self.getTransformation()

        #     # get parameters from other
        #     listOldParams = [param for _, param in oldParams.items()]
        #     listNewParams = [param for _, param in newParams.items()]

        #     def transAll(i, trans):
        #         stats_preprocessor = {}
        #         stats_preprocessor_reversed = {}

        #         for p, q in self.__stats_preprocessor.items():
        #             # change i-th value in p
        #             p = tuple([x if i != j else trans(p)[i] for j, x in enumerate(p)])

        #             # set new preprocessor stats
        #             stats_preprocessor[p] = q

        #             # check whether an error occured
        #             if stats_preprocessor_reversed.isContaining(q):
        #                 raise AttributeError('Internal error when the parameters are transformed to new parameter setting')
        #             stats_preprocessor_reversed[q] = p

        #         return stats_preprocessor, stats_preprocessor_reversed


        #     for i, (p1, p2) in enumerate(zip(listOldParams, listNewParams)):
        #         blub = False
        #         if p1.isUncertain() and not p2.isUncertain():
        #             # transform up (p1) -> dp (p2)
        #             s, sr = transAll(i, oldTrans.trans)
        #             blub = True
        #         elif not p1.isUncertain() and p2.isUncertain():
        #             # transform dp (p1) -> up (p2)
        #             s, sr = transAll(i, newTrans.inv_trans)
        #             blub = True

        #         if blub:
        #             print "-" * 60
        #             self.__stats_preprocessor = s
        #             self.__stats_preprocessor_reverse = sr

    # -----------------------------------------------------------------
    # UQSetting File formatter
    # -----------------------------------------------------------------
    def setMemento(self, memento):
        """
        Restores the state which is saved in the given memento
        @param memento: the memento object
        """
        self.fromJson(memento)

    def createMemento(self):
        """
        Creates a new memento to hold the current state
        @return: new memento
        """
        jsonString = self.toJson()
        jsonObject = json.loads(jsonString)
        return jsonObject

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the UQSetting object from the json object with its
        attributes.
        @param jsonObject: dictionary represeting a json object
        @return: the restored UQSetting object
        """
        setting = UQSetting()

        # restore samples
        key = '_UQSetting__stats_samples'
        if key in jsonObject:
            d = ju.parseKeyAsTuple(jsonObject[key])

            d0 = {}
            for key, value in d.items():
                d0[key] = Sample.fromJson(value)
            setting.setSamplesStats(d0)

        # restore preprocessor settings
        key = '_UQSetting__stats_preprocessor'
        if key in jsonObject:
            d1 = ju.parseKeyAsTuple(jsonObject[key])
            setting.setPreprocessorStats(d1)

        # restore postprocessor settings
        key = '_UQSetting__stats_preprocessor_reverse'
        if key in jsonObject:
            d2 = ju.parseKeyAsTuple(jsonObject[key])
            setting.setPreprocessorStatsReverse(d2)

        # restore simulation settings
        key = '_UQSetting__stats_simulation'
        if key in jsonObject:
            d3 = ju.parseKeyAsTuple(jsonObject[key])
            setting.setSimulationStats(d3)

        # restore postprocessor settings
        key = '_UQSetting__stats_postprocessor'
        if key in jsonObject:
            d = ju.parseKeyAsTuple(jsonObject[key])
            setting.setPostprocessorStats(d)

#             # --------------------------------------------------
#             # version 2
#             # --------------------------------------------------
#             dd = {}
#             # from pysgpp.extensions.datadriven.utils.GzipSerializer import GzipSerializer
#             broken_lmp, broken_reason, unbroken = [], [], []
#             print "-" * 60
#             for i, (q, _) in enumerate(d3.items()):
#                 isBroken = None
#                 # ----------------------------------------
#                 # read postlog
#                 # ----------------------------------------
#                 lmp = d3[q]['lmp']
#                 k = int(lmp.split('/')[-1])
#                 print "%i/%i (%i)" % (i + 1, len(d3), k)
#                 # dirname = os.path.dirname(lmp)
#                 dirname = '/home/franzefn/Promotion/Projekte/Peridynamik/results/results_fraunhofer/discretizations_uq3/sg'
#                 postlog = os.path.join(dirname, 'GCG-%i.postlog' % k)
#                 # post = os.path.join(dirname, 'sg_%i.post' % k)
#                 # ----------------------------------------
#                 # check postlog
#                 if not os.path.exists(postlog):
#                     isBroken = postlog
#                 else:
#                     # read file
#                     data = readDataARFF(postlog)['data'].array()
#                     if data.shape[0] != 1701:
#                         isBroken = '%i' % data.shape[0]
#                     else:
#                         # write them to the dictionary
#                         dd[q] = {}
#                         dd[q]['time'] = data[:, 0]
# #                             dd[q]['kinEng'] = data[:, 1]
# #                             dd[q]['potEng'] = data[:, 2]
# #                             dd[q]['te'] = data[:, 3]
# #                             dd[q]['dt'] = data[:, 4]
#                         dd[q]['gd'] = data[:, 5]
# #                     # read post
# #                     if not os.path.exists(post):
# #                         isBroken = post
# #                     else:
# #                         text = GzipSerializer().deserializeFromFile(post)
# #                         text = text.replace("'", '"')
# #                         jsonObject = json.loads(text)
# #                         dd[q]['force'] = jsonObject['force']
# #                         dd[q]['mean_damage'] = jsonObject['mean_damage']
# #                         dd[q]['var_damage'] = jsonObject['var_damage']
#
#                 # --------------------------------------------------
#                 # check if the parameter is broken
#                 # --------------------------------------------------
#                 if isBroken is not None:
#                     broken_lmp.append(lmp)
#                     broken_reason.append(isBroken)
#                     print "-" * 60
#                     print lmp
#                     print "-" * 60
#                 else:
#                     unbroken.append(lmp)
#
#             print "-" * 60
#             print len(broken_lmp), "+", len(unbroken), "=", \
#                 (len(broken_lmp) + len(unbroken))
#             print len(d0), len(d1), len(d2), len(d3), len(dd)
#
#             # write results to file
#             fd = open('broken_lmp.txt', 'w')
#             for lmp in broken_lmp:
#                 fd.write(lmp + "\n")
#             fd.close()
#
#             fd = open('broken_reason.txt', 'w')
#             for lmp in broken_reason:
#                 fd.write(lmp + "\n")
#             fd.close()
#
#             setting.setSamplesStats(d0)
#             setting.setPreprocessorStats(d1)
#             setting.setPreprocessorStatsReverse(d2)
#             setting.setSimulationStats(d3)
#             setting.setPostprocessorStats(dd)

            # --------------------------------------------------
            # version 1
            # --------------------------------------------------
            # print "-"*60
            # cnt_broken, cnt_unbroken = 0, 0
            # for p, item in d.items():
            #     if item.isContaining('_'):
            #         a = item['_'][0]['post']
            #         if not os.path.exists(a):
            #             continue #print 'file missing', a
            #         else:
            #             # read file
            #             fd = open(a, 'r')
            #             lines = fd.readlines()
            #             fd.close()

            #             if len(lines) > 0:
            #                 if len(lines) > 1:
            #                     lines = "".join(lines).replace('\n', ' ')
            #                 else:
            #                     lines = lines[0]

            #                 try:
            #                     l = lines.replace('nan', 'float("nan")')
            #                     l = lines.replace('inf', 'float("inf")')
            #                     dd[p] = eval(l)
            #                     if len(dd[p]['c_C1']) < 13:
            #                         print a
            #                         cnt_broken += 1
            #                     else:
            #                         cnt_unbroken += 1
            #                 except:
            #                     print a
            #                     cnt_broken += 1

            #     else:
            #         break
            # print "-"*60
            # print cnt_unbroken, "+", cnt_broken, "=", (cnt_broken + cnt_unbroken)
            # #if d[d.keys()[0]].isContaining('_'):
            # --------------------------------------------------
            # setting.setPostprocessorStats(dd)

        # assert len(setting.getPreprocessorStats()) ==  len(setting.getPreprocessorStatsReverse())
        # assert len(setting.getSimulationStats()) == len(setting.getPostprocessorStats())
        # assert len(setting.getPreprocessorStats()) == len(setting.getPostprocessorStats())

        # restore verbose setting
        key = '_UQSetting__verbose'
        if key in jsonObject:
            verbose = jsonObject[key] == 'true'
            setting.setVerbose(verbose)

        # restore last id
        key = 'lastid'
        if key in jsonObject:
            lastid = int(jsonObject[key])
            setting.setLastId(lastid)

        # restore specification settings
        spec = setting.getSpecification()

        # start time
        key = '_UQSpecification__t0'
        if key in jsonObject:
            spec.setStartTime(jsonObject[key])

        # end time
        key = '_UQSpecification__tn'
        if key in jsonObject:
            spec.setEndTime(jsonObject[key])

        # time step
        key = '_UQSpecification__dt'
        if key in jsonObject:
            spec.setTimeStep(jsonObject[key])

        # time step
        key = '_UQSpecification__reachesSteadyState'
        if key in jsonObject:
            spec.setReachesSteadyState(jsonObject[key])

        if len(setting.getPreprocessorStats()) != len(setting.getSimulationStats()):
            print 'Something went terribly wrong...'
            # import ipdb; ipdb.set_trace()

        # sanity check
        assert len(setting.getPreprocessorStats()) == len(setting.getSimulationStats())

        return setting

    def toJson(self):
        """
        @return: string representing the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in ('_UQSetting__stats_samples',
                         '_UQSetting__stats_preprocessor',
                         '_UQSetting__stats_preprocessor_reverse',
                         '_UQSetting__stats_simulation',
                         '_UQSetting__stats_postprocessor',
                         'lastid'):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(',\n')

        # print "j-------------------------------------------"
        # print "{" + s + "}"
        # print "j-------------------------------------------"

        return "{" + s + "}"

    def toDataMatrix(self, ts=[0], qoi='_',
                     dtype=SampleType.ACTIVEUNIT,
                     sampleType=UQSampleType.RAW,
                     *args, **kws):
        """
        Write simulation results to file
        @param ts: numeric time steps
        @param qoi: string quantity of interest
        """
        ans = {}
        res = self.getTimeDependentResults(ts, qoi,
                                           sampleType=sampleType,
                                           *args, **kws)
        for t, values in res.items():
            # if no number is specified here, just write all there is to file
            n = len(values)
            if n > 0:
                if sampleType == UQSampleType.RAW:
                    dim = values.iterkeys().next().getDim(dtype)
                else:
                    dim = len(values.iterkeys().next())

                data = DataMatrix(n, dim + 1)
                p = DataVector(dim + 1)

                # collect results
                for i, (key, res) in enumerate(values.items()):
                    if sampleType == UQSampleType.RAW:
                        key_value = key.getValue(dtype)
                    else:
                        key_value = key

                    for j, key_value_j in enumerate(key_value):
                        p[j] = float(key_value_j)
                    p[dim] = res
                    data.setRow(i, p)
                ans[t] = data
        return ans

    def __str__(self):
        return self.toJson()

    def defineSamples(self, params):
        """
        add samples object to the UQ setting if they are missing
        """
        for i, (p, _) in enumerate(self.__stats_preprocessor.items()):
            self.__stats_samples[p] = Sample(params, p, dtype=SampleType.EXPANDEDUNIT)

    def convert(self, params):
        # convert to new UQSetting
        # -> one needs to set the __stats_samples parameter
        if self.getSize("samples") < self.getSize("preprocessor"):
            self.defineSamples(params)

        for i, sample in enumerate(self.__stats_samples.values()):
            if params.getStochasticDim() != len(sample.getActiveUnit()):
                if i == 0:
                    warnings.warn("stochastic dimension changed -> applying it to the samples")
                # if the stochastic dimensionality has changed, let the samples know
                newSample = Sample(params, sample.getExpandedUnit(),
                                   SampleType.EXPANDEDUNIT)
                sample.init(newSample.getActiveUnit(),
                            newSample.getActiveProbabilistic(),
                            newSample.getExpandedUnit(),
                            newSample.getExpandedProbabilistic())
            # check if something changed in the accuracy of float
            p = tuple(sample.getExpandedUnit())
            if p not in self.__stats_samples:
                self.__stats_samples[p] = sample
                found = self.findEquivalent(sample, self.__stats_preprocessor)
                if not found:
                    raise AttributeError('can not find any results with respect to %s' % (p,))

    def cleanUp(self):
        """
        remove all the non-complete entries
        """
        for sample in self.__stats_samples.values():
            p = tuple(sample.getExpandedUnit())
            if p not in self.__stats_preprocessor:
                self.remove(sample)
            else:
                q = self.__stats_preprocessor[p]
                if q not in self.__stats_simulation or \
                        q not in self.__stats_postprocessor or \
                        len(self.__stats_postprocessor[q]) == 0:
                    self.remove(sample)

    def remove(self, sample):
        p = tuple(sample.getExpandedUnit())
        q = self.__stats_preprocessor[p]
        del self.__stats_samples[sample]
        del self.__stats_preprocessor[p]
        del self.__stats_preprocessor_reverse[q]
        del self.__stats_simulation[q]
        del self.__stats_postprocessor[q]

    def changeParamSetting(self, params):
        stats_samples = {}
        stats_preprocessor = {}
        stats_preprocessor_reverse = {}
        stats_simulation = {}
        stats_postprocessor = {}

        for p_old, sample_old in self.__stats_samples.items():
            # check if something changed in the accuracy of float
            p_new = np.array(p_old)
            hx = p_new[1]
            p_new[1] = p_new[2]
            p_new[2] = p_new[3]
            p_new[3] = hx
            sample_new = Sample(params, p_new, dtype=SampleType.EXPANDEDUNIT)
            assert np.all(sample_old.getActiveUnit() == sample_new.getActiveUnit())

            # load stats of old sample
            p_old = tuple(sample_old.getExpandedUnit())
            q_old = self.__stats_preprocessor[p_old]
            sim_old = self.__stats_simulation[q_old]
            post_old = self.__stats_postprocessor[q_old]

            p_new = tuple(sample_new.getExpandedUnit())
            q_new = np.array(q_old)
            hx = q_new[1]
            q_new[1] = q_new[2]
            q_new[2] = q_new[3]
            q_new[3] = hx
            q_new = tuple(q_new)

            # remove old one and add new one with old results
            stats_samples[p_new] = sample_new
            stats_preprocessor[p_new] = q_new
            stats_preprocessor_reverse[q_new] = p_new
            stats_simulation[q_new] = sim_old
            stats_postprocessor[q_new] = post_old

        self.__stats_samples = stats_samples
        self.__stats_preprocessor = stats_preprocessor
        self.__stats_preprocessor_reverse = stats_preprocessor_reverse
        self.__stats_simulation = stats_simulation
        self.__stats_postprocessor = stats_postprocessor

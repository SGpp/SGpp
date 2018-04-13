from multiprocessing import cpu_count
import inspect
import itertools
import math
import os
import pysgpp

from samplingresult import Samplingresult
import numpy as np
import remote_worker as remote


class UQSettingManager(object):
    """
    Interface to access the results stored in the UQSetting package
    """
    def __init__(self, uqsetting):
        self.uqsetting = uqsetting
        self.params = uqsetting.getParameters()
        self.setSampleGenerator()
        if not hasattr(self.uqsetting, 'lastid'):
            self.uqsetting.lastid = 0

        # parallel stuff
        if remote.do_distribution:
            self.parallelprocesses = sum(remote.hosts.itervalues())
            # make sure that the file is current
            self.uqsetting.writeToFile()
        else:
            # set to zero for no separate process
            self.parallelprocesses = cpu_count() / 2
        # set to expected number of samples, for parallelization
        self.expectedsamplecount = 0
        # incremented for each process
        self.__filesuffix = 0
        # running worker processes
        self.children = {}
        # result files written
        self.files = []

        if not self.uqsetting.getFilename():
            (_, filename, _, _, _, _) = inspect.getouterframes(inspect.currentframe())[1]
            self.uqsetting.setFilename(filename.replace(".py", "") +
                                       ".uqSetting.gz")
            self.uqsetting.writeToFile()

    def __next(self, number=1):
        """
        Allocate a new unique sample id. The id should be unique for the sample
        file, and increase monotonically.
        """
        idd = self.uqsetting.lastid + number
        self.uqsetting.lastid = idd
        return idd

    def getTags(self):
        """
        Find or create the tag list.
        """
        return self.uqsetting.getTagsStats()

    def sample(self, unitcubevalue, gen):
        """
        run uqsetting.run() and save the result
        @param unitcubevalue:the unit hyper cube sample to use
        @return: the sample like uqsetting.run()
        """
        if math.isnan(sum(unitcubevalue)):
            raise AttributeError("nan in unitcube!", unitcubevalue)
        t = gen.transform(gen.expandSample(unitcubevalue))
        if math.isnan(sum(t)):
            raise AttributeError("nan in unitcube!", t, unitcubevalue)
        self.uqsetting.run(t)
        # r['__unitcube_value'] = unitcubevalue
        return t

    def setExpectedSampleCount(self, n):
        self.expectedsamplecount = n

    def run_sampleList(self, sampleList, *tagList, **kwargs):
        """
        Performs run() for the given unit cube samples and tags them with the
        given tags. Parallelizes automatically, so use this if possible.
        Sampling will take place in the background, use waitForResults()
        at the end.
        """
        if not self.uqsetting.getFilename() in self.files:
            self.files.append(self.uqsetting.getFilename())

        if len(sampleList) == 0:
            return
        if self.parallelprocesses == 0:
            # just do the calculations
            self.do_sampleList(sampleList, tagList)
            return

        # split into smaller chunks
        jobsize = float(self.expectedsamplecount) / float(self.parallelprocesses)
        if jobsize < 1:
            jobsize = 1.0
        if len(sampleList) > jobsize:
            # round up
            njobs = math.ceil(len(sampleList) / jobsize)
            nsamples = int(math.ceil(len(sampleList) / njobs))
            njobs = int(njobs)
            print "jobconfig:", njobs, jobsize, nsamples, len(sampleList)

            # split into enough chunks to be below job size
            for i in range(0, njobs - 1):
                if (nsamples * (i + 1)) >= len(sampleList):
                    print i, "are enough"
                    break
                print "job", nsamples * i, \
                    len(sampleList[nsamples * i:nsamples * (i + 1)])
                self.run_sampleList(sampleList[nsamples * i:nsamples * (i + 1)],
                                    *tagList, starti=nsamples * i)
            # print "job", self.__starti, len(sampleList[nsamples*(njobs-1):])
            self.run_sampleList(sampleList[nsamples * (njobs - 1):],
                                *tagList, starti=nsamples * (njobs - 1))
            return

        if 'starti' in kwargs:
            starti = kwargs['starti']
        else:
            starti = self.__next(len(sampleList))

        # fork off one process per chunk
        if remote.do_distribution:
            host = remote.choose_host()
        else:
            host = ''

        filename = self.uqsetting.getFilename()\
                                 .replace(".gz", "." +
                                          str(self.__filesuffix) + ".gz")

        res = os.fork()

        if res == 0:
            # in the child, set a new file and do work

            if remote.do_distribution:
                remote.do_dist_run(host, self, filename, sampleList,
                                   tagList, starti)
                return  # unreachable

            self.uqsetting.setFilename(filename)
            self.uqsetting.setSimulationStats({})
            self.uqsetting.setPreprocessorStats({})
            self.uqsetting.setPreprocessorStatsReverse({})
            self.uqsetting.setPostprocessorStats({})
            self.uqsetting.setTagsStats({})
            self.do_sampleList(sampleList, tagList, starti)
            os._exit(0)  # die with no error

        elif res > 0:
            # in the parent, clean up and return
            self.children[res] = (starti, tagList, filename, host)
            # print "incremented suffix", map(lambda it: (it[0], it[1][0], it[1][2]), self.children.iteritems())
            self.__filesuffix = self.__filesuffix + 1
            return
        else:
            print "Error while forking, tags not processed:", tagList, starti

    def do_sampleList(self, sampleList, tagList, starti=0):
        for i, p in enumerate(sampleList):

            t = self.sample(p, self.gen)

            tags = self.uqsetting.getTagsStats()
            if t not in tags:
                tags[t] = []
            if tags[t] is None:
                tags[t] = []
            for tag in tagList:
                realtag = tag.copy()
                realtag['__index'] = starti + i
                tags[t].append(realtag)
        self.uqsetting.writeToFile()

    def waitForResults(self):
        """
        Wait for all forked workers to finish working.
        Then load all result files into one SamplingResult.
        """
        # this does not actually belong here, but here it will run in parallel.
        self.uqsetting.writeToFile()

        while len(self.children) > 0:
            pid, res = os.wait()
            if res == 0:
                self.files.append(self.children[pid][2])

                if self.children[pid][3] != '':
                    remote.free_host(self.children[pid][3])

                del self.children[pid]
                print "child finished, %d remaining" % len(self.children)
            else:
                starti, tags, f, host = self.children[pid]

                if host != '':
                    remote.free_host(host)

                print "WARNING: child %d crashed, File %s, Tags %s, starting from index %d on %s" % (pid, f, str(tags), starti, host)
                del self.children[pid]

    def loadResults(self):
        """
        Load Results from files, if necessary, and create a Samplingresult.
        """
        if self.files != []:
            if len(self.children) > 0:
                self.waitForResults()
            self.results = Samplingresult(self.files)
            self.files = []
            self.results.uq.setFilename(self.uqsetting.getFilename())
            self.uqsetting = self.results.uq
            self.uqsetting.writeToFile()
        else:
            self.results = Samplingresult(self.uqsetting)

    def getResults(self):
        """
        Return a SamplingResult with all the results.
        """
        if len(self.children) > 0:
            self.waitForResults()
            self.loadResults()
        else:
            self.loadResults()
        return self.results

    def run_distinct_paths(self, npaths, samples_per_path):
        # ---------------------------------------------------------------
        # Do monte-carlo UQ
        # ---------------------------------------------------------------

        # to get nice parallelisation
        self.setExpectedSampleCount(npaths * samples_per_path)

        # tag for anything that might be used for MC-analysis
        all_tag = {}
        all_tag['type'] = 'mc_all'
        all_tag['gen'] = self.gen.name

        genok = True

        for path in range(npaths):
            # restart the generator if necessary
            if not genok:
                self.gen.reset()

            # add a certain amount of random numbers
            pos = [self.gen.unitSample() for _ in range(samples_per_path)]
            genok = False

            # tag = self.createTag('distinct_mc_path', len=samples_per_path)
            chunk_tag = {}
            chunk_tag['type'] = 'mc_chunk'
            chunk_tag['id'] = path
            chunk_tag['gen'] = self.gen.name
            chunk_tag['len'] = samples_per_path
            # do the simulation
            self.run_sampleList(pos, chunk_tag, all_tag)

        self.waitForResults()

    def run_sensitivity(self, samples, maxDeg=1, samplingType='double'):
        # ---------------------------------------------------------------
        # Generate samples for ANOVA decomposition
        # Matrix A and B holding n samples, with k mixed matrix for
        # each parameter
        # samplingType: how to generate A and B.
        #   'restart':     restart the generator
        #   'double':      generate with twice the number of dimensions
        #   'interleave':  generate double number of samples and use even/odd
        #                  (note that this will actually create only half the number of sample points)
        # ---------------------------------------------------------------

        activeDim = self.gen.k

        maxDeg = min(activeDim, maxDeg)

        indexLists = []
        for deg in range(1, maxDeg + 1):
            indexLists = indexLists + [list(x) for x in itertools.combinations(range(activeDim), deg)]

        if samplingType == 'restart':
            samplesA = [self.gen.unitSample() for _ in range(samples)]
            self.gen.reset()
            samplesB = [self.gen.unitSample() for _ in range(samples)]
        elif samplingType == 'interleave':
            samples = samples / 2  # use no more samples as given by the user, to avoid problmes with sample generators.
            allSamples = [self.gen.unitSample() for _ in range(2 * samples)]
            samplesA = [allSamples[i] for i in range(0, 2 * samples, 2)]
            samplesB = [allSamples[i] for i in range(1, 2 * samples, 2)]
        elif samplingType == 'double':

            self.gen.k = 2 * activeDim
            self.gen.reset()
            allSamples = np.array([self.gen.unitSample()
                                   for _ in range(samples)])

            samplesA = allSamples[:, :activeDim]
            samplesB = allSamples[:, activeDim:]

            self.gen.k = activeDim
            self.gen.reset()

        else:
            print "ERROR: invalid samplingType!"

        # to get nice parallelisation
        self.setExpectedSampleCount(samples * (2 + len(indexLists)))

        # tag for anything that might be used for MC-analysis
        all_tag = {}
        all_tag['type'] = 'mc_all'
        all_tag['gen'] = self.gen.name
        all_tag['samplingType'] = samplingType

        sensitivity_tag = {}
        sensitivity_tag['type'] = 'mc_sensitivity'
        sensitivity_tag['gen'] = self.gen.name
        sensitivity_tag['len'] = samples
        sensitivity_tag['samplingType'] = samplingType

        # Matrix A generieren
        sensitivity_tag['part'] = 'A'
        self.run_sampleList(samplesA, sensitivity_tag, all_tag)

        # Matrix B generieren
        # samplesB = [self.gen.unitSample() for _ in range(samples)]
        sensitivity_tag['part'] = 'B'
        self.run_sampleList(samplesB, sensitivity_tag, all_tag)

        for indices in indexLists:
            sensitivity_tag['part'] = str(indices)
            print str(indices)
            self.run_sampleList(self.__mixMatrices(samplesA, samplesB, indices),
                                sensitivity_tag)

        self.waitForResults()

    def setSampleGenerator(self, generator=pysgpp.NaiveSampleGenerator):
        self.gen = Sampler(self.params, generator)

    def __mixMatrices(self, mat1, mat2, indices):
        result = []
        for i in range(len(mat1)):
            sample = list(mat1[i])
            for k in indices:
                sample[k] = mat2[i][k]
            result.append(tuple(sample))
        return result

    def __powerset(self, lst):
        f = lambda result, x: result + [subset + [x] for subset in result]
        return reduce(f, lst, [[]])


class Sampler:

    def __init__(self, parameterset, generator=pysgpp.NaiveSampleGenerator):
        self.n = parameterset.getDim()
        self.relevant = [k for k, i in parameterset.items() if i.isActive()]
        self.k = len(self.relevant)

        if callable(generator):
            self.gen = generator(self.k)
            self.gengen = generator
        else:
            self.gen = generator
            self.gengen = None

        self.params = parameterset
        self.name = self.gen.__class__.__name__

    def reset(self):
        if self.gengen:
            self.gen = self.gengen(self.k)
        else:
            print "WARNING: cannot reset generator"

    def unitSample(self):
        dv = pysgpp.DataVector(self.k)
        self.gen.getSample(dv)
        return [dv[i] for i in range(self.k)]

    def expandSample(self, unitsample):
        r = [-1] * self.n
        i = 0
        for i, j in zip(self.relevant, range(self.k)):
            r[i] = unitsample[j]
        return tuple(r)

    def transform(self, fullunitsample):
        t = [None] * self.n
        for i, p in self.params.items():
            if not p.isActive():
                t[i] = p.getSample()
            else:
                t[i] = p.getDistribution().ppf(fullunitsample[i])
        return self.params.trans(tuple(t))

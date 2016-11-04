class BuildUQSetting(object):

    def defineUQSetting(self, uqSettingFile):
        raise NotImplementedError

    def defineParameters(self, numDims, setting=None):
        raise NotImplementedError
    
    def __init__(self,
                 name,
                 param_setting,
                 collocation_settings):
        self.radix = 'test_parabola'

        # change working directory if necessary
        self.oldcwd = os.getcwd()
        if self.oldcwd.find('parabola') == -1:
            os.chdir(os.path.join(self.oldcwd, 'parabola'))

        # set up object variables
        self.numDims = 2
        self.param_setting = "uniform"

        # create folder for results
        self.pathResults = os.path.join('results/%s' % self.param_setting)


        self.defineParameters(self.numDims, self.param_setting)

        # available labels
        labels = ['sg', 'scc', 'fg', 'ref']
        filenames = ["%s.%s.%s.uqSetting.gz" % (self.radix, self.param_setting, label)
                     for label in labels]
        self.uqSettingsFilenames = dict(zip(labels, filenames))

        # define UQSettings
        self.uqSettings = {}
        for label, filename in self.uqSettingsFilenames.items():
            print "Read %s" % filename,
            self.uqSettings[label] = self.defineUQSetting(filename)
            print self.uqSettings[label].getSize(), \
                self.uqSettings[label].getAvailableQoI()
            self.uqSettings[label].convert(self.params)

        # compute reference values
        if 'ref' in self.uqSettings:
            self.computeReferenceValues(self.uqSettings['ref'])

    def computeReferenceValues(self, uqSetting, n=10000):
        # ----------------------------------------------------------
        # analytic reference values
        # ----------------------------------------------------------
        g = uqSetting.getSimulation()
        numDims = self.params.getStochasticDim()
        U = self.params.getIndependentJointDistribution()
        computeWithMC = False
        if self.param_setting == "uniform":
            print "computing analytic results"
            self.E_ana = (2. / 3.) ** numDims, 0.0
            if numDims == 1:
                self.V_ana = 4. / 45., 0.0
            elif numDims == 2:
                self.V_ana = 176. / 2025., 0.0
            elif numDims == 3:
                self.V_ana = 60416. / 820125., 0.0
            elif numDims == 4:
                self.V_ana = 1705984. / 36905625., 0.0
            else:
                computeWithMC = True
        else:
            if numDims == 1:
                print "computing analytic results 1d"
                self.E_ana = quad(lambda x: g([x]) * U.pdf([x]), 0, 1)
                self.V_ana = quad(lambda x: (g([x]) - self.E_ana[0]) ** 2 * U.pdf([x]), 0, 1)
            elif numDims == 2:
                print "computing analytic results 2d"
                self.E_ana = dblquad(lambda x, y: g([x, y]) * U.pdf([x, y]),
                                    0, 1, lambda x: 0, lambda x: 1)
                self.V_ana = dblquad(lambda x, y: (g([x, y]) - self.E_ana[0]) ** 2 * U.pdf([x, y]),
                                    0, 1, lambda x: 0, lambda x: 1)
            else:
                computeWithMC = True

        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        print "computing monte carlo reference values"
        n -= uqSetting.getSize()
        if n > 0:
            mcSampler = MCSampler.withLatinHypercubeSampleGenerator(self.params, n)
            samples = mcSampler.nextSamples(n)
            uqSetting.runSamples(samples)
            uqSetting.writeToFile()

        # ----------------------------------------------------------
        # monte carlo reference values
        # ----------------------------------------------------------
        res = uqSetting.getResults()
        analysis = MCAnalysis(self.params, res)

        if computeWithMC:
            print "computing analytic results > 2d"
            self.E_ana = analysis.mean()
            self.V_ana = analysis.var()

        self.refSize = len(res)

        # ----------------------------------------------
        # write reference values to file
        # ----------------------------------------------
        analysis.writeMoments("results/%s/%s.mc" % (self.param_setting, self.param_setting))

        # write reference values to file
        stats = {'data': [[self.E_ana[0]],
                           [self.E_ana[1]],
                           [self.V_ana[0]],
                           [self.V_ana[1]]],
                 'names': ["mean", "meanError", "var", "varError"],
                 'filename': "results/%s/%s.ref.moments.arff" % (self.param_setting, self.param_setting)}
        writeDataARFF(stats)

        print "-" * 60
        print "E(f) = %.14f, %g" % self.E_ana
        print "V(f) = %.14f, %g" % self.V_ana
        print "-" * 60

    def buildSetting(self, label, level, deg,
                     clenshaw_curtis=False,
                     modified=False,
                     nsamples=1000,
                     isFull=False,
                     epsilon=1e-15,
                     adaptive=None,
                     knowledgeFilename=None):
        builder = ASGCUQManagerBuilder()

        builder.withParameters(self.params)\
               .useUQSetting(self.uqSettings[label])\
               .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                      KnowledgeTypes.SQUARED])\
               .useInterpolation()

        if 'ref' in self.uqSettings and len(self.uqSettings['ref']) > 0:
            builder.withTestSet(self.uqSettings['ref'])\
                   .learnWithTest()\

        if knowledgeFilename is not None:
            builder.withKnowledge(knowledgeFilename)

        samplerSpec = builder.defineSampler()
        gridSpec = samplerSpec.withGrid()
        gridSpec.withLevel(level)
        if deg > 1:
            gridSpec.withPolynomialBase(deg)
        if isFull:
            gridSpec.isFull()
        if clenshaw_curtis:
            gridSpec.isClenshawCurtis()
        if modified:
            gridSpec.withModifiedBasis()

        if adaptive is not None:
            # specify the refinement
            samplerSpec.withRefinement()\
                       .withAdaptThreshold(epsilon)\
                       .withAdaptPoints(2)\
                       .withBalancing()

            refinement = samplerSpec.refineMostPromisingNodes()
            if adaptive == "simple":
                refinement.withSurplusRanking()
            elif adaptive == "exp":
                refinement.withExpectationValueOptimizationRanking()
            elif adaptive == "var":
                refinement.withVarianceOptimizationRanking()
            elif adaptive == "squared":
                refinement.withSquaredSurplusRanking()

            refinement.createAllChildrenOnRefinement()

            samplerSpec.withStopPolicy().withGridSizeLimit(nsamples)

        uqManager = builder.andGetResult()

        # update the stats, which are not stored with the knowledge
        if knowledgeFilename is not None:
            uqManager.recomputeStats()

        return uqManager


    def runSampler(self, uqManager, label, alabel, blabel):
        uqSetting = self.uqSettings[label]
        # ----------------------------------------------
        # prepare folders
        pathResults = os.path.join(self.pathResults, alabel, blabel)

        if not os.path.exists(pathResults):
            os.makedirs(pathResults)

        for newdir in [os.path.join(pathResults, 'checkpoints'),
                       os.path.join(pathResults, 'grids'),
                       os.path.join(pathResults, 'samples')]:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        # write the setting to file
        uqManager.uqSetting.writeToFile()

    def defineASGCAnalysis(self, uqManager):
        builder = ASGCAnalysisBuilder()
        builder.withUQManager(uqManager)\
               .withAnalyticEstimationStrategy()
        return builder

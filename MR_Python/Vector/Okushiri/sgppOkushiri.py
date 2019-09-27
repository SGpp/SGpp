import execnet
import numpy as np
import pysgpp
import pickle


def call_python_version(Version, Module, Function, ArgumentList):
    # Use execnet to call ANUGA, which is still only python2 from a python3 environment
    gw = execnet.makegateway("popen//python=python%s" % Version)
    channel = gw.remote_exec("""
        from %s import %s as the_function
        channel.send(the_function(*channel.receive()))
    """ % (Module, Function))
    channel.send(ArgumentList)
    result = channel.receive()
    gw.exit()
    return result


class okushiriStorage():
    # hashes evaluations of the okushiri model, to speed up SG++ runtimes
    def __init__(self, dim, numTimeSteps, gridResolution):
        self.dim = dim
        self.numTimeSteps = numTimeSteps
        self.gridResolution = gridResolution
        np.savetxt(
            '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/numTimeSteps.txt', [numTimeSteps])
        np.savetxt(
            '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/gridsize.txt', [gridResolution])
        self.precalcValuesFileName = '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/precalculatedValues/sg_precalculations{}D{}T{}R.pkl'.format(
            dim, numTimeSteps, gridResolution)
        try:
            with open(self.precalcValuesFileName, 'rb') as f:
                self.precalculatedValues = pickle.load(f)
        except:
            self.precalculatedValues = {}
        self.numNew = 0

    def cleanUp(self):
        with open(self.precalcValuesFileName, "wb") as f:
            pickle.dump(self.precalculatedValues, f)
        print("\ncalculated {} new Okushiri evaluations".format(self.numNew))
        if self.numNew > 0:
            print("saved them to {}".format(self.precalcValuesFileName))

    def eval(self, x):
        parameters = np.zeros(self.dim)
        for i in range(self.dim):
            parameters[i] = x[i]
        # lists are not allowed as keys, but tuples are
        key = tuple(parameters)
        if key in self.precalculatedValues:
            y = self.precalculatedValues[key]
        else:
            print("sgppOkushiri: processing {}".format(parameters))
            # x = np.array([0.5, 0.5, 0.875, 0.5])
            np.savetxt(
                '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/x.txt', parameters)

            # reset time
            np.savetxt(
                '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/t.txt', [-1])

            # run solver
            # This must currently be called from /home/rehmemk/git/SGpp/MR_Python/Vector
            call_python_version("2.7", "Okushiri.okushiri", "run", [])
            y = np.loadtxt(
                '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/data/y.txt')
            #print("Ran okushiri benchmark with python 2.7 in {}s".format(dt))

            # plt.plot(range(len(y)), y)
            # plt.show()

            self.precalculatedValues[key] = y
            self.numNew += 1
            print("sgppOkushiri: Done ({})".format(self.numNew))
            # if self.numNew % 50 == 0:
            #     self.cleanUp()
        return y


class okushiri():
    # The Okushiri Benchmark
    def __init__(self, dim, numTimeSteps=451):
        self.dim = dim
        self.out = numTimeSteps
        self.gridResolution = 16
        self.okushiriStorage = okushiriStorage(
            dim, numTimeSteps, self.gridResolution)
        self.pdfs = pysgpp.DistributionsVector()
        for _ in range(dim):
            self.pdfs.push_back(pysgpp.DistributionUniform(0.0, 1.0))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "okushiri{}D{}T{}R".format(self.getDim(), self.getOut(), self.gridResolution)

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.out

    def eval(self, x):
        y = self.okushiriStorage.eval(x)
        return y

    def cleanUp(self):
        self.okushiriStorage.cleanUp()

    def evalJacobian(self, x):
        print("Jacobian unknown")
        return 0

    def getDistributions(self):
        return self.pdfs

    def getMeans(self):
        print("Means unknown")
        means = np.zeros(self.getOut())
        return means

    def getVars(self):
        print("Vars unknown")
        variances = np.zeros(self.getOut())
        return variances

    def getIntegrals(self):
        print("Integrals unknown")
        integrals = np.zeros(self.getOut())
        return integrals

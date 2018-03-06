class ASGCStatistics(object):

    def __init__(self):
        # statistics per knowledge type
        self.trainMSE = {}
        self.trainL2Norm = {}
        self.trainCount = {}

        self.testMSE = {}
        self.testL2Norm = {}
        self.testL1Norm = {}
        self.testMaxError = {}
        self.testCount = {}

        self.numberPoints = {}
        self.level = {}
    
    def __updateDicts(self, dtype, t):
        if dtype not in self.trainMSE:
            self.trainMSE[dtype] = {}
            self.trainL2Norm[dtype] = {}
            self.trainCount[dtype] = {}

            self.testMSE[dtype] = {}
            self.testL2Norm[dtype] = {}
            self.testL1Norm[dtype] = {}
            self.testMaxError[dtype] = {}
            self.testCount[dtype] = {}

            self.numberPoints[dtype] = []
            self.level[dtype] = []
            
        if t not in self.trainMSE[dtype]:
            self.trainMSE[dtype][t] = []
            self.trainL2Norm[dtype][t] = []
            self.trainCount[dtype][t] = []

            self.testMSE[dtype][t] = []
            self.testL2Norm[dtype][t] = []
            self.testL1Norm[dtype][t] = []
            self.testMaxError[dtype][t] = []
            self.testCount[dtype][t] = []
    
    def updateResults(self, dtype, t, learner):
        # update the dictionaries
        self.__updateDicts(dtype, t)

        # evaluate MSE of training data set -> should be zero
        self.trainMSE[dtype][t].append(learner.getMSE("train"))
        self.trainL2Norm[dtype][t].append(learner.getL2NormError("train"))
        self.trainCount[dtype][t].append(learner.getSize("train"))

        # store interpolation quality
        if learner.getSize("test") > 0:
            # L2 error and MSE
            self.testMSE[dtype][t].append(learner.getMSE("test"))
            self.testL2Norm[dtype][t].append(learner.getL2NormError("test"))
            self.testL1Norm[dtype][t].append(learner.getL1NormError("test"))
            self.testMaxError[dtype][t].append(learner.getMaxError("test"))
            self.testCount[dtype][t].append(learner.getSize("test"))

        # update other stats
        gs = learner.grid.getStorage()
        self.level[dtype].append(gs.getMaxLevel())
        self.numberPoints[dtype].append(gs.getSize())

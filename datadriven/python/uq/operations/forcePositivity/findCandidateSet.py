class CandidateSet(object):


    def __init__(self):
        self.candidates_history = {}
        self.candidates = []
        self.costs = 0
        self.iteration = 0
        self.verbose = True


    def findGridPointsWithNegativeCoefficient(self, gps, alpha):
        negativeGridPoints = {}
        for i, gp in gps.items():
            if alpha[i] < 0.0:
                negativeGridPoints[i] = gp
        return negativeGridPoints


    def separateGridPoints(self, gps, alpha):
        negativeGridPoints = {}
        positiveGridPoints = {}
        for i, gp in gps.items():
            if alpha[i] < 0.0:
                negativeGridPoints[i] = gp
            else:
                positiveGridPoints[i] = gp

        return negativeGridPoints, positiveGridPoints


    def removeAlreadyExistingGridPoints(self, grid, intersections):
        gs = grid.getStorage()
        return [gp for gp in intersections if not gs.has_key(gp)]

    
    def findCandidates(self):
        raise NotImplementedError


    def hasMoreCandidates(self, grid, alpha, addedGridPoints):
        self.findCandidates(grid, alpha, addedGridPoints)
        self.candidates_history[self.iteration] = self.candidates
        self.iteration += 1
        return len(self.candidates) > 0
    
    
    def nextCandidateSet(self):
        ans = self.candidates, self.costs
        self.candidates = []
        self.costs = 0
        return ans


# ----------------------------------------------------------------------------

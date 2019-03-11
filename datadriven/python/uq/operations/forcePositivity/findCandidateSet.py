class CandidateSet(object):


    def __init__(self):
        self.candidates_history = {}
        self.candidates = []
        self.costs = 0
        self.iteration = 0
        self.verbose = True

    def removeAlreadyExistingGridPoints(self, grid, intersections):
        gs = grid.getStorage()
        return [gp for gp in intersections if not gs.isContaining(gp)]


    def findCandidates(self, grid, alpha, addedGridPoints):
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

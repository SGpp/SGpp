from spatiallyAdaptiveBase import *


class SpatiallyAdaptiveSingleDimensions(SpatiallyAdaptivBase):
    def __init__(self, a, b, grid=None):
        SpatiallyAdaptivBase.__init__(self, a, b, grid)

    def coarsenGrid(self, area, levelvec):
        pass
        '''
        #a = self.a[0]
        #b = self.b[0]
        containedInGrid = True #the areas might be very fine -> it might not be contained in a coarse grid
        startOfGridRange = True #this indicates that the area starts at a gridpoint on the grid
        startArea = np.zeros(self.dim)
        endArea = np.zeros(self.dim)
        for d in range(self.dim):
            startArea[d] = self.refinement[d][area[d]][0]
            endArea[d] = self.refinement[d][area[d]][1]
            width = (endArea[d] - startArea[d])
            if (self.b[d]-self.a[d])/width > 2**levelvec[d]: #this means area is too fine
                containedInGrid = False
                #check if area is at start of points in the grid -> extend area to grid resolution
                #this avoids that we extend areas more than once
                widthCoarse = (self.b[d]-self.a[d])/float(2**levelvec[d])
                if int(startArea[d]/widthCoarse) != startArea[d] / widthCoarse:
                    startOfGridRange = False
                else:
                    endArea[d] = startArea[d] + widthCoarse
        if(not containedInGrid and not startOfGridRange): #this area is not considered in this grid
            return None, startArea, endArea
        #calculate the number of points in this area
        levelCoarse = np.zeros(self.dim,dtype=int)
        for d in range(self.dim):
            coarseningFactor = self.refinement[d][area[d]][2]
            levelCoarse[d] = max(int(levelvec[d]  +  math.log(((endArea[d] - startArea[d]) / (self.b[d]-self.a[d])*1.0), 2)   - coarseningFactor), 0)
        return levelCoarse, startArea, endArea
        '''

    # returns the points of a single component grid with refinement
    def get_points_component_grid(self, levelvec, numSubDiagonal):
        pass
        '''
        array2 = []
        areas = self.getAreas()
        #iterate through all of the subareas to generate list of points
        for area in areas:
            levelCoarse, start, end = self.coarsenGrid(area,levelvec)
            if(levelCoarse == None):
                continue
            #generate the points of the area; same concept as before with indices -> generate tuples of size dim with all point combinations
            self.grid.setCurrentArea(start,end,levelCoarse)
            points = self.grid.getPoints()
            array2.extend(points)
        return array2
        '''

    # this method draws the 1D refinement of each dimension individually
    def drawRefinement(self, filename=None):  # update with meta container
        plt.rcParams.update({'font.size': 32})
        refinement = self.refinement
        dim = len(refinement)
        fig, ax = plt.subplots(ncols=1, nrows=dim, figsize=(20, 10))
        for d in range(dim):
            starts = [refinementObject.start for refinementObject in refinement.refinementContainers[d].get_objects()]
            ends = [refinementObject.end for refinementObject in refinement.refinementContainers[d].get_objects()]
            for i in range(len(starts)):
                ax[d].add_patch(
                    patches.Rectangle(
                        (starts[i], -0.1),
                        ends[i] - starts[i],
                        0.2,
                        fill=False  # remove background
                    )
                )
            xValues = starts + ends
            yValues = np.zeros(len(xValues))
            ax[d].plot(xValues, yValues, 'bo', markersize=10, color="black")
            ax[d].set_xlim([self.a[d], self.b[d]])
            ax[d].set_ylim([-0.1, 0.1])
            ax[d].set_yticks([])
        if filename is not None:
            plt.savefig(filename, bbox_inches='tight')
        plt.show()
        return fig

    # evaluate the integral of f in a specific area with numPoints many points using the specified integrator set in the grid
    # We also interpolate the function to the finest width to calculate the error of the combination in each
    def evaluateArea(self, f, area, levelvec):
        pass
        # redo with hierarchization
        '''
        levelCoarse, start, end = self.coarsenGrid(area,levelvec)
        if(levelCoarse == None):
            #return default value to indicate that we do not compute this area
            return -2**30, [], 0
        startArea = np.zeros(self.dim)
        endArea = np.zeros(self.dim)
        for d in range(self.dim):
            startArea[d] = self.refinement[d][area[d]][0]
            endArea[d] = self.refinement[d][area[d]][1]
        width = endArea - startArea
        partialIntegral = self.grid.integrate(f,levelCoarse,start,end) #calculate regular integral
        partialIntegrals = [(partialIntegral, start, end)]
        numGridsOfRefinement = np.ones(self.dim)
        #calculate the number of grids in each dimension that we need to interpolate for the error estimator
        if(not(np.equal(np.array(end)-np.array(start),np.array(width)).all())):
            numGridsOfRefinement = [int((end[d] - start[d]) / self.finestWidth )for d in range(self.dim)]
        #calculate the integral of the interpolated partial grids
        if(np.prod(numGridsOfRefinement) > 1):
            self.grid.setCurrentArea(start,end,levelvec)
            partialIntegrals = integrateVariableStartNumpyArbitraryDimAndInterpolate(f,self.grid.levelToNumPoints(levelCoarse),start,end,numGridsOfRefinement)
        return partialIntegral, partialIntegrals, np.prod(self.grid.levelToNumPoints(levelCoarse))

        '''

    # calculate the position of the partial integral in the variance array
    def getPosition(self, partialIntegralInfo, dim):
        refinement_array = self.refinement
        start_partial_integral = partialIntegralInfo[1]
        end_partial_integral = partialIntegralInfo[2]
        positions = []
        for d in range(dim):
            position_dim = []
            for i in range(len(refinement_array[d])):
                r = refinement_array[d][i]
                start_ref = r[0]
                end_ref = r[1]
                if start_partial_integral[d] >= start_ref and end_partial_integral[d] <= end_ref:
                    position_dim.append(i)
            positions.append(position_dim)
        return positions

    def initializeRefinement(self):
        # toDo self, start, end, dim, coarseningLevel = 0
        RefinementObjectSingleDimension
        initial_points = []
        for d in range(self.dim):
            initial_points.append(np.linspace(self.a[d], self.b[d], 2 ** 1 + 1))
        self.refinement = MetaRefinementContainer([RefinementContainer
                                                   ([RefinementObjectSingleDimension(initial_points[d][i],
                                                                                     initial_points[d][i + 1], self.dim,
                                                                                     self.lmax[d] - 1) for i in
                                                     range(2 ** 1)], self.dim, self.errorEstimator) for d in
                                                   range(self.dim)])
        # self.finestWidth = (initial_points[0][1]-initial_points[0][0])/(self.b[0] - self.a[0])

    def getAreas(self):
        # get a list of lists which contains range(refinements[d]) for each dimension d where the refinements[d] are the number of subintervals in this dimension
        indices = [list(range(len(refineDim))) for refineDim in self.refinement.get_new_objects()]
        # this command creates tuples of size dim of all combinations of indices (e.g. dim = 2 indices = ([0,1],[0,1,2,3]) -> areas = [(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)] )
        return list(zip(*[g.ravel() for g in np.meshgrid(*indices)]))

    def getNewAreas(self):
        return self.get_areas()

    def prepareRefinement(self):
        pass
        # this is done in meta container
        '''
        self.popArray = [[] for d in range(self.dim)]
        self.newScheme = False
        '''

    def doRefinement(self, area, position):
        lmaxChange = self.refinement.refine(position)
        if lmaxChange is not None:
            self.lmax = [self.lmax[d] + lmaxChange[d] for d in range(self.dim)]
            print("New scheme")
            self.scheme = self.combischeme.getCombiScheme(self.lmin[0], self.lmax[0], self.dim)
            return True
        return False

    def refinementPostprocessing(self):
        '''
        if self.newScheme:
            #restore scheme to initial state
            initialPoints = np.linspace(self.a[0],self.b[0],2**1+1)
            self.refinement = [[(initialPoints[i],initialPoints[i+1],self.lmax[d]-1) for i in range(2**1)] for d in range(self.dim)]
            self.finestWidth = initialPoints[1] - initialPoints[0]
        else:
            #remove outdated refinement values
            for d in range(self.dim):
                for position in reversed(sorted(list(set(self.popArray[d])))):
                    self.refinement[d].pop(position)
        '''
        self.refinement.apply_remove()
        self.refinement.reinit_new_objects()
        # self.combiintegral = 0.0

    '''
    def getErrors(self,integralarrayComplete, errorOperator, f):
        #errorArray = np.zeros(len(self.getAreas()))
        offsetAreas = np.ones(self.dim, dtype=int)
        for d in range(self.dim - 1):
            offsetAreas[d+1] = len(self.refinement[d]) * offsetAreas[d]
        for i in range(len(integralarrayComplete)):
            for j in range(len(integralarrayComplete[i][0])):
                position = self.getPosition(integralarrayComplete[i][0][j],self.dim)
                error = errorOperator.calcError(f,integralarrayComplete[i][0][j])
                for d in range(self.dim):
                    assert(len(position[d]) == 1)
                p = np.inner(np.ravel(position), offsetAreas)
                #errorArray[p] += error * integralarrayComplete[i][1]
        #self.errorArray = list(errorArray)
    '''

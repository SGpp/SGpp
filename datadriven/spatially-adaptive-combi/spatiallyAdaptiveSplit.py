from spatiallyAdaptiveBase import *


class SpatiallyAdaptiveFixedScheme(SpatiallyAdaptivBase):

    # returns the points of a single component grid with refinement
    def get_points_component_grid(self, levelvec, numSubDiagonal):
        dim = len(levelvec)
        points_array = []
        for area in self.refinement.get_objects():
            start = area.start
            end = area.end
            levelInterval = np.zeros(dim, dtype=int)
            for d in range(dim):
                levelInterval[d] = int(levelvec[d] - self.lmin[d])
            self.grid.setCurrentArea(start, end, levelInterval)
            points = self.grid.getPoints()
            points_array.extend(points)
        # print points_array
        return points_array

    def get_points_and_weights_component_grid(self, levelvec, numSubDiagonal):
        dim = len(levelvec)
        points_array = []
        weights_array = []
        for area in self.refinement.get_objects():
            start = area.start
            end = area.end
            levelInterval = np.zeros(dim, dtype=int)
            for d in range(dim):
                levelInterval[d] = int(levelvec[d] - self.lmin[d])
            self.grid.setCurrentArea(start, end, levelInterval)
            points, weights = self.grid.get_points_and_weights()
            points_array.extend(points)
            weights_array.extend(weights)
        # print array2
        return points_array, weights_array

    # draw a visual representation of refinement tree
    def draw_refinement(self, filename=None):
        plt.rcParams.update({'font.size': 32})
        dim = len(self.refinement[0][0])
        if dim > 2:
            print("Refinement can only be printed in 2D")
            return
        fig = plt.figure(figsize=(20, 20))
        ax2 = fig.add_subplot(111, aspect='equal')
        # print refinement
        for i in self.refinement.get_objects():
            startx = i.start[0]
            starty = i.start[1]
            endx = i.end[0]
            endy = i.end[1]
            ax2.add_patch(
                patches.Rectangle(
                    (startx, starty),
                    endx - startx,
                    endy - starty,
                    fill=False  # remove background
                )
            )
            # content = str(i[5])
            xcontent = startx + (endx - startx) / 2.0
            ycontent = starty + (endy - starty) / 2.0
            # print content, xcontent, ycontent
            # ax.annotate(content, (xcontent,ycontent))
        if filename is not None:
            plt.savefig(filename, bbox_inches='tight')
        plt.show()
        return fig

    def initialize_refinement(self):
        parent_integral = self.grid.integrate(self.f, np.zeros(self.dim, dtype=int), self.a, self.b)
        parent = RefinementObjectExtendSplit(np.array(self.a), np.array(self.b), self.grid,
                                             self.numberOfRefinementsBeforeExtend, None, 0,
                                             0)
        parent.set_integral(parent_integral)
        new_refinement_objects = parent.split_area_arbitrary_dim()
        self.refinement = RefinementContainer(new_refinement_objects, self.dim, self.errorEstimator)
        if self.errorEstimator is None:
            self.errorEstimator = ErrorCalculatorExtendSplit()

    def evaluate_area(self, f, area, levelvec):
        levelEval = np.zeros(self.dim, dtype=int)
        for d in range(self.dim):
            levelEval[d] = int(levelvec[d] - self.lmin[d])
        return self.grid.integrate(f, levelEval, area.start, area.end), None, np.prod(
            self.grid.levelToNumPoints(levelEval))

    def do_refinement(self, area, position):
        self.refinement.refine(position)
        return False

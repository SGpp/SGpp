from spatiallyAdaptiveBase import *


class SpatiallyAdaptiveExtendScheme(SpatiallyAdaptivBase):
    def __init__(self, a, b, number_of_refinements_before_extend=1, grid=None, no_initial_splitting=False,
                 version=0, dim_adaptive=False, automatic_extend_split=False):
        # there are three different version that coarsen grids slightly different
        # version 0 coarsen as much as possible while extending and adding only new points in regions where it is supposed to
        # version 1 coarsens less and also adds moderately many points in non refined regions which might result in a more balanced configuration
        # version 2 coarsen fewest and adds a bit more points in non refinded regions but very similar to version 1
        assert 2 >= version >= 0
        self.version = version
        SpatiallyAdaptivBase.__init__(self, a, b, grid)
        self.noInitialSplitting = no_initial_splitting
        self.numberOfRefinementsBeforeExtend = number_of_refinements_before_extend
        self.refinements_for_recalculate = 100
        self.dim_adaptive = dim_adaptive
        self.automatic_extend_split = automatic_extend_split

    # draw a visual representation of refinement tree
    def draw_refinement(self, filename=None):
        plt.rcParams.update({'font.size': 32})
        dim = self.dim
        if dim > 2:
            print("Refinement can only be printed in 2D")
            return
        fig = plt.figure(figsize=(20, 20))
        ax2 = fig.add_subplot(111, aspect='equal')
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
        if filename is not None:
            plt.savefig(filename, bbox_inches='tight')
        plt.show()
        return fig

    # returns the points of a single component grid with refinement
    def get_points_component_grid(self, levelvec, numSubDiagonal):
        assert (numSubDiagonal < self.dim)
        points_array = []
        for area in self.refinement.get_objects():
            start = area.start
            end = area.end
            level_interval, is_null = self.coarsen_grid(levelvec, area, numSubDiagonal)
            self.grid.setCurrentArea(start, end, level_interval)
            points = self.grid.getPoints()
            points_array.extend(points)
        return points_array

    def get_points_and_weights_component_grid(self, levelvec, numSubDiagonal):
        assert (numSubDiagonal < self.dim)
        points_array = []
        weights_array = []
        for area in self.refinement.get_objects():
            start = area.start
            end = area.end
            level_interval, is_null = self.coarsen_grid(levelvec, area, numSubDiagonal)
            self.grid.setCurrentArea(start, end, level_interval)
            points, weights = self.grid.get_points_and_weights()
            points_array.extend(points)
            weights_array.extend(weights)
        return points_array, weights_array

    # returns the points of a single component grid with refinement
    def get_points_component_grid_not_null(self, levelvec, numSubDiagonal):
        assert (numSubDiagonal < self.dim)
        array2 = []
        for area in self.refinement.get_objects():
            start = area.start
            end = area.end
            level_interval, is_null = self.coarsen_grid(levelvec, area, numSubDiagonal)
            if not is_null:
                self.grid.setCurrentArea(start, end, level_interval)
                points = self.grid.getPoints()
                array2.extend(points)
                # print("considered", levelvec, level_interval, area.start, area.end, area.coarseningValue)
            # else:
            # print("not considered", levelvec, level_interval, area.start, area.end, area.coarseningValue)
        return array2

    # optimized adaptive refinement refine multiple cells in close range around max variance (here set to 10%)
    def coarsen_grid(self, levelvector, area, num_sub_diagonal, print_point=None):
        start = area.start
        end = area.end
        coarsening = area.coarseningValue
        temp = list(levelvector)
        coarsening_save = coarsening
        area_is_null = False
        if self.version == 0:

            maxLevel = max(temp)
            temp2 = list(reversed(sorted(list(temp))))
            if temp2[0] - temp2[1] < coarsening:
                while coarsening > 0:
                    maxLevel = max(temp)
                    if maxLevel == self.lmin[0]:  # we assume here that lmin is equal everywhere
                        break
                    for d in range(self.dim):
                        if temp[d] == maxLevel:
                            temp[d] -= 1
                            coarsening -= 1
                            break
                area_is_null = True
            else:
                for d in range(self.dim):
                    if temp[d] == maxLevel:
                        temp[d] -= coarsening
                        break
                if area.is_already_calculated(tuple(temp), tuple(levelvector)):
                    area_is_null = True
                else:
                    area.add_level(tuple(temp), tuple(levelvector))
        else:
            while coarsening > 0:
                maxLevel = max(temp)
                if maxLevel == self.lmin[0]:  # we assume here that lmin is equal everywhere
                    break
                occurences_of_max = 0
                for i in temp:
                    if i == maxLevel:
                        occurences_of_max += 1
                is_top_diag = num_sub_diagonal == 0
                if self.version == 1:
                    no_forward_problem = coarsening_save >= self.lmax[0] + self.dim - 1 - maxLevel - (
                            self.dim - 2) - maxLevel + 1
                    do_coarsen = no_forward_problem and coarsening >= occurences_of_max - is_top_diag
                else:
                    no_forward_problem = coarsening_save >= self.lmax[0] + self.dim - 1 - maxLevel - (
                            self.dim - 2) - maxLevel + 2
                    do_coarsen = no_forward_problem and coarsening >= occurences_of_max
                if do_coarsen:
                    for d in range(self.dim):
                        if temp[d] == maxLevel:
                            temp[d] -= 1
                            coarsening -= 1
                else:
                    break
        level_coarse = [temp[d] - self.lmin[d] + int(self.noInitialSplitting) for d in range(len(temp))]
        if print_point is not None:
            if all([start[d] <= print_point[d] and end[d] >= print_point[d] for d in range(self.dim)]):
                print("Level: ", levelvector, "Coarsened level:", level_coarse, coarsening_save, start, end)
        return level_coarse, area_is_null

    def initialize_refinement(self):
        if self.dim_adaptive:
            self.combischeme.init_adaptive_combi_scheme(self.lmax, self.lmin)
        if self.noInitialSplitting:
            new_refinement_object = RefinementObjectExtendSplit(np.array(self.a), np.array(self.b), self.grid,
                                                                self.numberOfRefinementsBeforeExtend, 0, 0, automatic_extend_split=self.automatic_extend_split)
            self.refinement = RefinementContainer([new_refinement_object], self.dim, self.errorEstimator)
        else:
            parent_integral = self.grid.integrate(self.f, np.zeros(self.dim, dtype=int), self.a, self.b)
            parent = RefinementObjectExtendSplit(np.array(self.a), np.array(self.b), self.grid,
                                                                 self.numberOfRefinementsBeforeExtend, None, 0,
                                                                 0, automatic_extend_split=self.automatic_extend_split)
            parent.set_integral(parent_integral)
            new_refinement_objects = parent.split_area_arbitrary_dim()
            self.refinement = RefinementContainer(new_refinement_objects, self.dim, self.errorEstimator)
        if self.errorEstimator is None:
            self.errorEstimator = ErrorCalculatorExtendSplit()

    def evaluate_area(self, f, area, levelvec):
        num_sub_diagonal = (self.lmax[0] + self.dim - 1) - np.sum(levelvec)
        level_for_evaluation, is_null = self.coarsen_grid(levelvec, area, num_sub_diagonal)
        if is_null:
            return 0, None, 0
        else:
            return self.grid.integrate(f, level_for_evaluation, area.start, area.end), None, np.prod(
                self.grid.levelToNumPoints(level_for_evaluation))

    def do_refinement(self, area, position):
        if self.automatic_extend_split:
            if area.extend_parent_integral is None:
                lmax_reduced = [self.lmax[d] - 1 for d in range(self.dim)]
                area.extend_parent_integral = 0.0
                scheme = self.combischeme.getCombiScheme(self.lmin[0], lmax_reduced[0], self.dim, False)
                for ss in scheme:
                    area_integral, partial_integrals, evaluations = self.evaluate_area(self.f, area, ss[0])
                    area.extend_parent_integral += area_integral * ss[1]
            if area.split_parent_integral is None:
                area.parent.coarseningValue = area.coarseningValue
                area.parent.levelvec_dict = {}
                area.split_parent_integral = 0
                for ss in self.scheme:
                    area_integral, partial_integrals, evaluations = self.evaluate_area(self.f, area.parent, ss[0])
                    area.split_parent_integral += area_integral * ss[1]

        lmax_change = self.refinement.refine(position)
        if lmax_change != None:
            self.lmax = [self.lmax[d] + lmax_change[d] for d in range(self.dim)]
            print("New scheme")
            self.scheme = self.combischeme.getCombiScheme(self.lmin[0], self.lmax[0], self.dim)
            return True
        return False

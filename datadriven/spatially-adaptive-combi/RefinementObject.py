import math
import numpy as np
import abc, logging
from combiScheme import CombiScheme


# This class defines the general template of a Refinement Object that can be stored in the refinement container
class RefinementObject(object):
    def __init__(self, error_estimator):
        self.errorEstimator = error_estimator
        self.integral = None

    # set the local integral for area associated with RefinementObject
    def set_integral(self, integral):
        self.integral = integral

    # refine this object and return newly created objects
    @abc.abstractmethod
    def refine(self):
        pass

    # set the local error associated with RefinementObject
    def set_error(self, error):
        self.error = error


# This is the special class for the RefinementObject defined in the split extend scheme
class RefinementObjectExtendSplit(RefinementObject):
    def __init__(self, start, end, grid, number_of_refinements_before_extend, parent_integral, coarseningValue=0,
                 needExtendScheme=0, punish_depth=False, extend_parent_integral=None,
                 split_parent_integral=None, automatic_extend_split=False, parent=None, factor=1,
                 depth=0):
        # start of subarea
        self.start = start
        # end of subarea
        self.end = end
        self.dim = len(start)
        # indicates how often area needs to be coarsened according to extend scheme
        self.coarseningValue = coarseningValue
        # indicates how many splits where already performed for this region
        self.needExtendScheme = needExtendScheme
        # indicates after how many splits only extends are performed
        self.numberOfRefinementsBeforeExtend = number_of_refinements_before_extend
        self.evaluations = 0
        self.integral = None
        self.parent_integral = parent_integral
        # dictionary that maps all coarsened levelvectors to there uncoarsened ones
        # the can only be one uncoarsened levelvector for each coarsened one all other areas are set to 0
        self.levelvec_dict = {}
        self.punish_depth = punish_depth
        self.grid = grid
        self.error_extend = None
        self.error_split = None
        self.automatic_extend_split = automatic_extend_split
        self.extend_parent_integral =  extend_parent_integral
        self.split_parent_integral = split_parent_integral
        self.parent = parent
        self.children = []
        self.factor = factor
        self.depth = depth

    # this routine decides if we split or extend the RefinementObject
    def refine(self):
        coarsening_level = self.coarseningValue
        if self.automatic_extend_split:
            #print(self.start, self.end, "Parent", self.parent_integral, "Integral", self.integral, "Extend",
            #      self.extend_parent_integral, "Split", self.split_parent_integral)
            if self.error_extend is None:
                #print(self.extend_parent_integral, self.integral)
                self.error_extend = abs(self.extend_parent_integral - self.integral) / self.dim
            if self.error_split is None:
                sum_siblings = 0.0
                i = 0
                for child in self.parent.children:
                    if child.integral is not None:
                        sum_siblings += child.integral
                        i += 1
                #print(i)
                assert i == 2 ** self.dim  # we always have 2**dim children
                self.error_split = abs(self.split_parent_integral - sum_siblings) / (
                            2 ** self.dim * 2 ** (self.depth ** 2))  # 2**self.dim)
                # self.error_split = abs(self.split_parent_integral/2**self.dim - self.integral)/ 2**(self.depth)#math.sqrt(2**self.dim)
            #print("Extend error", self.error_extend, "Split error", self.error_split)
        if (self.automatic_extend_split and self.error_extend > self.error_split) or (
                not self.automatic_extend_split and
                self.needExtendScheme >= self.numberOfRefinementsBeforeExtend):  # add new component grids to scheme and refine only target area

            if self.coarseningValue == 0:
                coarseningValue = 0
            else:
                coarseningValue = coarsening_level - 1
            # in case we have refined complete scheme (i.e. coarensingLevel was 0)
            # we have to increase level everywhere else
            if coarsening_level == 0:
                # increase lmax by dim
                lmaxIncrease = [1 for d in range(self.dim)]
                newRefinementObject = RefinementObjectExtendSplit(self.start, self.end, self.grid,
                                                                  self.numberOfRefinementsBeforeExtend, self.integral,
                                                                  coarseningValue, self.needExtendScheme,
                                                                  extend_parent_integral=self.integral,
                                                                  automatic_extend_split=self.automatic_extend_split,
                                                                  parent=self.parent, factor=self.dim,
                                                                  depth=self.depth + 0.5 * math.sqrt(self.dim))
                self.children.append(newRefinementObject)
                return [newRefinementObject], lmaxIncrease, 1
            else:
                # add to integralArray
                newRefinementObject = RefinementObjectExtendSplit(self.start, self.end, self.grid,
                                                                  self.numberOfRefinementsBeforeExtend, self.integral,
                                                                  coarseningValue, self.needExtendScheme,
                                                                  extend_parent_integral=self.integral,
                                                                  automatic_extend_split=self.automatic_extend_split,
                                                                  parent=self.parent, factor=self.dim,
                                                                  depth=self.depth + 0.5 * math.sqrt(self.dim))
                self.children.append(newRefinementObject)
                return [newRefinementObject], None, None
        elif (self.automatic_extend_split and self.error_extend <= self.error_split) or (
                not self.automatic_extend_split and self.needExtendScheme >= 0):  # split the array
            # add to integralArray
            self.needExtendScheme += 1
            #print("Splitting", self.start, self.end)
            newRefinementObjects = self.split_area_arbitrary_dim()
            return newRefinementObjects, None, None
        else:
            print("Error!!!! Invalid value")
            assert False

    # in case lmax was changed the coarsening value of other RefinementObjects need to be increased
    def update(self, update_info):
        self.coarseningValue += update_info
        self.levelvec_dict = {}

    def add_level(self, levelvec_coarsened, levelvec):
        # print("adding", levelvec_coarsened, levelvec)
        self.levelvec_dict[levelvec_coarsened] = levelvec

    def is_already_calculated(self, levelvec_coarsened, levelvec):
        if levelvec_coarsened not in self.levelvec_dict:
            return False
        else:
            # print(levelvec_coarsened, levelvec, self.levelvec_dict[levelvec_coarsened])
            return self.levelvec_dict[levelvec_coarsened] != levelvec

    # splits the current area into 2**dim smaller ones and returns them
    def split_area_arbitrary_dim(self):
        dim = self.dim
        num_sub_areas = 2 ** dim
        start = self.start
        end = self.end
        midpoint = [self.grid.get_mid_point(start[d], end[d], d) for d in range(self.dim)]
        sub_area_array = []
        for i in range(num_sub_areas):
            start_sub_area = np.zeros(dim)
            end_sub_area = np.zeros(dim)
            rest = i
            for d in reversed(list(range(dim))):
                start_sub_area[d] = start[d] if rest < 2 ** d else midpoint[d]
                end_sub_area[d] = midpoint[d] if rest < 2 ** d else end[d]
                rest = rest % 2 ** d
            new_refinement_object = RefinementObjectExtendSplit(start_sub_area, end_sub_area, self.grid,
                                                                self.numberOfRefinementsBeforeExtend,
                                                                self.integral / 2 ** self.dim,
                                                                self.coarseningValue, self.needExtendScheme,
                                                                split_parent_integral=self.integral,
                                                                parent=self,
                                                                automatic_extend_split=self.automatic_extend_split,
                                                                factor=1, depth=self.depth + 0.25 * math.sqrt(self.dim))
            self.children.append(new_refinement_object)
            sub_area_array.append(new_refinement_object)
        return sub_area_array

    # set the local error associated with RefinementObject
    def set_error(self, error):
        if (self.punish_depth):
            error = error * np.prod(np.array(self.end) - np.array(self.start)) * 2 ** self.coarseningValue
        self.error = error


# This is the special class for the RefinementObject defined in the split extend scheme
class RefinementObjectCell(RefinementObject):
    cell_dict = {}
    cells_for_level = []
    punish_depth = False

    def __init__(self, start, end, levelvec, a, b, lmin, father=None):
        self.a = a
        self.b = b
        self.lmin = lmin
        # start of subarea
        self.start = start
        # end of subarea
        self.end = end
        RefinementObjectCell.cell_dict[self.get_key()] = self
        self.dim = len(start)
        self.levelvec = np.array(levelvec, dtype=int)
        # print("levelvec", self.levelvec)
        self.level = sum(levelvec) - self.dim + 1
        # self.father = father
        self.children = []
        # self.descendants = set()
        self.error = None
        self.active = True
        self.parents = []
        self.sub_integrals = []
        self.integral = None
        for d in range(self.dim):
            parent = RefinementObjectCell.parent_cell_arbitrary_dim(d, self.levelvec, self.start, self.end, self.a,
                                                                    self.b, self.lmin)
            if parent is not None:
                self.parents.append(parent)
                if parent != father:
                    # print(parent, RefinementObjectCell.cell_dict.items(), self.get_key(), father, self.levelvec)
                    parent_object = RefinementObjectCell.cell_dict[parent]
                    parent_object.add_child(self)

    def add_child(self, child):
        self.children.append(child)

    def get_key(self):
        return tuple((tuple(self.start), tuple(self.end)))

    def isActive(self):
        return self.active

    def contains(self, point):
        contained = True
        for d in range(self.dim):
            if point[d] < self.start[d] or point[d] > self.end[d]:
                contained = False
                break
        return contained

    def is_corner(self, point):
        is_corner = True
        for d in range(self.dim):
            if point[d] != self.start[d] and point[d] != self.end[d]:
                is_corner = False
                break
        return is_corner

    def refine(self):
        assert self.active
        self.active = False
        new_objects = []
        for d in range(self.dim):
            levelvec_copy = list(self.levelvec)
            levelvec_copy[d] += 1
            possible_candidates_d = RefinementObjectCell.children_cell_arbitrary_dim(d, self.start, self.end, self.dim)
            for candidate in possible_candidates_d:
                if candidate in RefinementObjectCell.cell_dict:
                    continue
                # print("candidate", candidate)
                # key = candidate.get_key()
                can_be_refined = True
                for parent in RefinementObjectCell.get_parents(levelvec_copy, candidate[0], candidate[1], self.a,
                                                               self.b, self.dim, self.lmin):
                    if parent not in RefinementObjectCell.cell_dict or RefinementObjectCell.cell_dict[
                        parent].isActive():
                        can_be_refined = False
                        break
                if can_be_refined:
                    new_objects.append(
                        RefinementObjectCell(candidate[0], candidate[1], list(levelvec_copy), self.a, self.b, self.lmin,
                                             self.get_key()))

        self.children.extend(new_objects)
        # print("New refined objects", [object.get_key() for object in new_objects])
        return new_objects, None, None

    # splits the current cell into the 2 children in the dimension d and returns the children
    def split_cell_arbitrary_dim(self, d):
        childs = RefinementObjectCell.children_cell_arbitrary_dim(d, self.start, self.end, self.dim)
        sub_area_array = []
        levelvec = list(self.levelvec)
        levelvec[d] += 1
        for child in childs:
            if child not in RefinementObjectCell.cell_dict:
                new_refinement_object = RefinementObjectCell(child[0], child[1], list(levelvec), self.a, self.b,
                                                             self.lmin, self.get_key())
                # RefinementObjectCell.cells_for_level[self.level+1].append(new_refinement_object)
                sub_area_array.append(new_refinement_object)
        return sub_area_array

    # splits the current cell into the 2 children in the dimension d and returns the children
    @staticmethod
    def children_cell_arbitrary_dim(d, start, end, dim):
        spacing = np.zeros(dim)
        spacing[d] = 0.5 * (end[d] - start[d])
        start_sub_area = np.array(start)
        end_sub_area = np.array(end)
        child1 = tuple((tuple(start_sub_area + spacing), tuple(end_sub_area)))
        child2 = tuple((tuple(start_sub_area), tuple(end_sub_area - spacing)))
        return [child1, child2]

    @staticmethod
    def get_parents(levelvec, start, end, a, b, dim, lmin):
        parents = []
        for d in range(dim):
            parent = RefinementObjectCell.parent_cell_arbitrary_dim(d, levelvec, start, end, a, b, lmin)
            if parent is not None:
                parents.append(parent)
        return parents

    @staticmethod
    def parent_cell_arbitrary_dim(d, levelvec, start, end, a, b, lmin):
        levelvec_parent = list(levelvec)
        if levelvec_parent[d] <= lmin[d]:
            return None
        levelvec_parent[d] = levelvec_parent[d] - 1
        parent_start = np.array(start)
        parent_end = np.array(end)
        index_of_start = start[d] * 2 ** levelvec[d]
        if index_of_start % 2 == 1:  # start needs to be changed
            parent_start[d] = parent_end[d] - (b[d] - a[d]) / 2 ** levelvec_parent[d]
        else:  # end needs to be changed
            parent_end[d] = parent_start[d] + (b[d] - a[d]) / 2 ** levelvec_parent[d]

        return tuple((tuple(parent_start), tuple(parent_end)))

    def get_points(self):
        return set(zip(*[g.ravel() for g in np.meshgrid(*[[self.start[d], self.end[d]] for d in range(self.dim)])]))

    # sets the error of a cell only if it is refinable; otherwise we do not need to set it here
    def set_error(self, error):
        # only define an error if active cell
        if self.active:
            self.error = error
            if self.punish_depth:
                self.error *= np.prod(np.array(self.end) - np.array(self.start))
        else:
            self.error = 0
        # print("Error of refine object:", self.get_key(), "is:", self.error)

        self.sub_integrals = []


# This is the special class for the RefinementObject defined in the split extend scheme
class RefinementObjectSingleDimension(RefinementObject):
    def __init__(self, start, end, dim, coarsening_level=0):
        # start of subarea
        self.start = start
        # end of subarea
        self.end = end
        self.dim = dim
        # indicates how often area needs to be coarsened according to extend scheme
        self.coarsening_level = coarsening_level
        self.evaluations = 0

    def refine(self):
        # coarseningLevel = self.refinement[dimValue][area[dimValue]][2]
        # if coarseningLevel == 0 and allDimsMaxRefined == False: #not used currently
        #    continue
        # positionDim = area[dimValue]
        # if positionDim in self.popArray[dimValue]: #do not refine one interval twice
        #    continue
        # print("refine:", positionDim , "in dim:", dimValue, "coarsening Level:", coarseningLevel)
        # splitAreaInfo = self.refinement[dimValue][positionDim]
        # lowerBorder = splitAreaInfo[0]
        # upperBorder = splitAreaInfo[1]
        # self.popArray[dimValue].append(positionDim)
        lmax_increase = None
        update = None
        coarsening_value = 0
        # if we are already at maximum refinement coarsening_value stays at 0
        if self.coarsening_level == 0:
            coarsening_value = 0
        else:  # otherwise decrease coarsening level
            coarsening_value = self.coarsening_level - 1
        # in case we have refined complete scheme (i.e. coarensingLevel was 0) we have to increase level everywhere else
        if (self.coarsening_level == 0):  # extend scheme if we are at maximum refinement
            # increase coarsening level in all other dimensions
            for d in range(self.dim):
                for i in range(len(self.refinement[d])):
                    icoarsen = self.refinement[d][i][2] + 1
                    self.refinement[d][i] = (self.refinement[d][i][0], self.refinement[d][i][1], icoarsen)
            # increase lmax by 1
            lmax_increase = [1 for d in range(self.dim)]
            update = 1
            # print("New scheme")
            # self.scheme = getCombiScheme(self.lmin[0],self.lmax[0],self.dim)
            # self.newScheme = True
        # add new refined interval to refinement array (it has half of the width)
        new_width = (self.end - self.start) / 2.0
        new_objects = []
        new_objects.append(
            RefinementObjectSingleDimension(self.start, self.start + new_width, self.dim, coarsening_value))
        new_objects.append(
            RefinementObjectSingleDimension(self.start + new_width, self.end, self.dim, coarsening_value))
        # self.finestWidth = min(new_width,self.finestWidth)
        return new_objects, lmax_increase, update

    # in case lmax was changed the coarsening value of other RefinementObjects need to be increased
    def update(self, update_info):
        self.coarsening_level += update_info

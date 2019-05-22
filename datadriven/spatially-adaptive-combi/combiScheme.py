import numpy as np
import math


class CombiScheme:
    def __init__(self, dim):
        self.initialized_adaptive = False
        self.active_index_set = set()
        self.old_index_set = set()
        self.dim = dim

    def init_adaptive_combi_scheme(self, lmax, lmin):
        assert lmax >= lmin
        self.lmin = lmin
        self.lmax = lmax
        self.initialized_adaptive = True
        self.active_index_set = CombiScheme.init_active_index_set(lmax, lmin, self.dim)
        self.old_index_set = CombiScheme.init_old_index_set(lmax, lmin, self.dim)
        self.lmax_adaptive = lmax

    def extendable_level(self, levelvec):
        assert self.initialized_adaptive
        counter = 0
        extendable_dim = 0
        for d in range(self.dim):
            if levelvec[d] > 1:
                counter += 1
                extendable_dim = d
        return counter == 1, extendable_dim

    def is_refinable(self, levelvec):
        assert self.initialized_adaptive
        return tuple(levelvec) in self.active_index_set

    def update_adaptive_combi(self, levelvec):
        assert self.initialized_adaptive
        if not self.is_refinable(levelvec):
            return
        refined_dims = []
        # remove this levelvec from active_index_set and add to old_index_set
        self.active_index_set.remove(tuple(levelvec))
        self.old_index_set.add(tuple(levelvec))
        for d in range(self.dim):
            if self.__refine_scheme(d, levelvec):
                refined_dims.append(d)
        return refined_dims

    def __refine_scheme(self, d, levelvec):
        assert self.initialized_adaptive
        # print(CombiScheme.old_index_set, CombiScheme.active_index_set, levelvec, CombiScheme.lmin)
        levelvec = list(levelvec)
        levelvec[d] += 1
        for dim in range(self.dim):
            levelvec_copy = list(levelvec)
            levelvec_copy[dim] = levelvec[dim] - 1
            if tuple(levelvec_copy) not in self.old_index_set and not levelvec_copy[dim] < self.lmin:
                return False
        self.active_index_set.add(tuple(levelvec))
        self.lmax_adaptive = max(self.lmax_adaptive, levelvec[d])
        return True

    @staticmethod
    def init_active_index_set(lmax, lmin, dim):
        grids = CombiScheme.getGrids(dim, lmax - lmin + 1)
        grids = [tuple([l + (lmin - 1) for l in g]) for g in grids]
        #print(grids)
        return set(grids)

    @staticmethod
    def init_old_index_set(lmax, lmin, dim):
        grid_array = []
        for q in range(1, min(dim, lmax - lmin + 1)):
            grids = CombiScheme.getGrids(dim, lmax - lmin + 1 - q)
            grids = [tuple([l + (lmin - 1) for l in g]) for g in grids]
            grid_array.extend(grids)
        #print(grid_array)
        return set(grid_array)

    def get_index_set(self):
        return self.old_index_set | self.active_index_set

    def getCombiScheme(self, lmin, lmax, dim, do_print=False):
        grid_array = []
        if not self.initialized_adaptive:  # use default scheme
            for q in range(min(dim, lmax - lmin + 1)):
                coefficient = (-1) ** q * math.factorial(dim - 1) / (math.factorial(q) * math.factorial(dim - 1 - q))
                grids = CombiScheme.getGrids(dim, lmax - lmin + 1 - q)
                grid_array.extend(
                    [(np.array(g, dtype=int) + np.ones(dim, dtype=int) * (lmin - 1), coefficient) for g in grids])
            for i in range(len(grid_array)):
                if do_print:
                    print(i, list(grid_array[i][0]), grid_array[i][1])
        else:  # use adaptive schem
            assert self.initialized_adaptive
            grid_array = self.get_coefficients_to_index_set(self.active_index_set | self.old_index_set)
            # print(grid_dict.items())
            for i in range(len(grid_array)):
                if do_print:
                    print(i, list(grid_array[i][0]), grid_array[i][1])
        return grid_array

    def get_coefficients_to_index_set(self, index_set):
        grid_array = []
        grid_dict = {}
        for grid_levelvec in index_set:
            stencils = []
            for d in range(self.dim):
                if grid_levelvec[d] <= self.lmin:
                    stencils.append([0])
                else:
                    stencils.append([0, -1])
            stencil_elements = list(zip(*[g.ravel() for g in np.meshgrid(*stencils)]))
            for s in stencil_elements:
                levelvec = tuple(map(lambda x, y: x + y, grid_levelvec, s))  # adding tuples
                update_coefficient = -(abs((sum(s))) % 2) + (abs(((sum(s)) - 1)) % 2)
                if levelvec in grid_dict:
                    grid_dict[levelvec] += update_coefficient
                else:
                    grid_dict[levelvec] = update_coefficient

        for levelvec, coefficient in grid_dict.items():
            if coefficient != 0 and sum(levelvec) > self.lmax + (self.dim - 1) * self.lmin - self.dim:
                grid_array.append((levelvec, coefficient))
        return grid_array

    def is_old_index(self, levelvec):
        return tuple(levelvec) in self.old_index_set

    @staticmethod
    def getGrids(dim_left, values_left):
        if dim_left == 1:
            return [[values_left]]
        grids = []
        for index in range(values_left):
            levelvector = [index + 1]
            grids.extend([levelvector + g for g in CombiScheme.getGrids(dim_left - 1, values_left - index)])
        return grids

    def in_index_set(self, levelvec):
        return tuple(levelvec) in self.active_index_set or tuple(levelvec) in self.old_index_set

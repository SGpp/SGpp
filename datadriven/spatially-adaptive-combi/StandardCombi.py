import matplotlib.pyplot as plt
import abc
from Function import FunctionShift
import logging
from combiScheme import *

# T his class implements the standard combination technique
class StandardCombi(object):
    # initialization
    # a = lower bound of integral; b = upper bound of integral
    # grid = specified grid (e.g. Trapezoidal);
    def __init__(self, a, b, grid=None):
        self.log = logging.getLogger(__name__)
        self.dim = len(a)
        self.a = a
        self.b = b
        self.grid = grid
        self.combischeme = CombiScheme(self.dim)
        assert (len(a) == len(b))

    def set_combi_parameters(self, minv, maxv):
        # compute minimum and target level vector
        self.lmin = [minv for i in range(self.dim)]
        self.lmax = [maxv for i in range(self.dim)]
        # get combi scheme
        self.scheme = self.combischeme.getCombiScheme(minv, maxv, self.dim)

    # standard combination scheme for quadrature
    # lmin = minimum level; lmax = target level
    # f = function to integrate; dim=dimension of problem
    def perform_combi(self, minv, maxv, f, reference_solution=None):
        start = self.a
        end = self.b
        self.f = f
        self.f.reset_dictionary()
        self.set_combi_parameters(minv, maxv)
        combiintegral = 0
        for ss in self.scheme:
            integral = self.grid.integrate(self.f, ss[0], start, end) * ss[1]
            combiintegral += integral
        real_integral = reference_solution
        print("CombiSolution", combiintegral)
        if reference_solution is not None:
            print("Analytic Solution", real_integral)
            print("Difference", abs(combiintegral - real_integral))
            return self.scheme, abs(combiintegral - real_integral), combiintegral
        else:
            return self.scheme, None, combiintegral

    def get_num_points_component_grid(self, levelvector, doNaive, num_sub_diagonal):
        return np.prod(self.grid.levelToNumPoints(levelvector))

    # calculate the total number of points used in the complete combination scheme
    def get_total_num_points(self, doNaive = False,
                                           distinct_function_evals=True):  # we assume here that all lmax entries are equal
        if distinct_function_evals:
            return self.f.get_f_dict_size()
        numpoints = 0
        for ss in self.scheme:
            num_sub_diagonal = (self.lmax[0] + self.dim - 1) - np.sum(ss[0])
            pointsgrid = self.get_num_points_component_grid(ss[0], doNaive, num_sub_diagonal)
            if distinct_function_evals and self.grid.isNested():
                numpoints += pointsgrid * int(ss[1])
            else:
                numpoints += pointsgrid
        # print(numpoints)
        return numpoints

    # prints every single component grid of the combination and orders them according to levels
    def print_resulting_combi_scheme(self, filename=None):
        plt.rcParams.update({'font.size': 22})
        scheme = self.scheme
        lmin = self.lmin
        lmax = [self.combischeme.lmax_adaptive]*self.dim if hasattr(self.combischeme, 'lmax_adaptive') else self.lmax
        dim = self.dim
        if dim != 2:
            print("Cannot print combischeme of dimension > 2")
            return None
        fig, ax = plt.subplots(ncols=lmax[0] - lmin[0] + 1, nrows=lmax[1] - lmin[1] + 1, figsize=(20, 20))
        # get points of each component grid and plot them individually
        if lmax == lmin:
            ax.xaxis.set_ticks_position('none')
            ax.yaxis.set_ticks_position('none')
            ax.set_xlim([self.a[0] - 0.005, self.b[0] + 0.005])
            ax.set_ylim([self.a[1] - 0.005, self.b[1] + 0.005])
            num_sub_diagonal = (self.lmax[0] + dim - 1) - np.sum(lmax)
            points = self.get_points_component_grid(lmax, num_sub_diagonal)
            x_array = [p[0] for p in points]
            y_array = [p[1] for p in points]
            ax.plot(x_array, y_array, 'o', markersize=6, color="black")
        else:
            for i in range(lmax[0] - lmin[0] + 1):
                for j in range(lmax[1] - lmin[1] + 1):
                    ax[i, j].xaxis.set_ticks_position('none')
                    ax[i, j].yaxis.set_ticks_position('none')
                    ax[i, j].set_xlim([self.a[0] - 0.005, self.b[0] + 0.005])
                    ax[i, j].set_ylim([self.a[1] - 0.005, self.b[1] + 0.005])
            for ss in scheme:
                num_sub_diagonal = (self.lmax[0] + dim - 1) - np.sum(ss[0])
                points = self.get_points_component_grid(ss[0], num_sub_diagonal)
                x_array = [p[0] for p in points]
                y_array = [p[1] for p in points]
                ax[lmax[1] - lmin[1] - (ss[0][1] - lmin[1]), (ss[0][0] - lmin[0])].plot(x_array, y_array, 'o', markersize=6,
                                                                                        color="black")
        if filename is not None:
            plt.savefig(filename, bbox_inches='tight')
        plt.show()
        return fig

    # prints the sparse grid which results from the combination
    def print_resulting_sparsegrid(self, filename=None):
        plt.rcParams.update({'font.size': 32})
        scheme = self.scheme
        dim = self.dim
        if dim != 2:
            print("Cannot print sparse grid of dimension > 2")
            return None
        fig, ax = plt.subplots(figsize=(20, 20))
        ax.set_xlim([self.a[0] - 0.01, self.b[0] + 0.01])
        ax.set_ylim([self.a[1] - 0.01, self.b[1] + 0.01])
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        # ax.axis('off')
        # get points of each component grid and plot them in one plot
        for ss in scheme:
            numSubDiagonal = (self.lmax[0] + dim - 1) - np.sum(ss[0])
            points = self.get_points_component_grid(ss[0], numSubDiagonal)
            xArray = [p[0] for p in points]
            yArray = [p[1] for p in points]
            plt.plot(xArray, yArray, 'o', markersize=10, color="black")
        if filename is not None:
            plt.savefig(filename, bbox_inches='tight')
        plt.show()
        return fig

    # check if combischeme is right
    def check_combi_scheme(self):
        if not self.grid.isNested():
            return
        dim = self.dim
        dictionary = {}
        for ss in self.scheme:
            num_sub_diagonal = (self.lmax[0] + dim - 1) - np.sum(ss[0])
            # print num_sub_diagonal , ii ,ss
            points = self.get_points_component_grid_not_null(ss[0], num_sub_diagonal)
            points = set(points)
            for p in points:
                if p in dictionary:
                    dictionary[p] += ss[1]
                else:
                    dictionary[p] = ss[1]
        # print(dictionary.items())
        for key, value in dictionary.items():
            # print(key, value)
            if value != 1:
                print(dictionary)
                print("Failed for:", key, " with value: ", value)
                for area in self.refinement.get_objects():
                    print("area dict", area.levelvec_dict)
                '''
                for area in self.refinement.getObjects():
                    print("new area:",area)
                    for ss in self.scheme:
                        num_sub_diagonal = (self.lmax[0] + dim - 1) - np.sum(ss[0])
                        self.coarsenGrid(ss[0],area, num_sub_diagonal,key)
                #print(self.refinement)
                #print(dictionary.items())
                '''
            assert (value == 1)

    def get_points_component_grid_not_null(self, levelvec, numSubDiagonal):
        return self.get_points_component_grid(levelvec, numSubDiagonal)

    def get_points_component_grid(self, levelvec, numSubDiagonal):
        self.grid.setCurrentArea(self.a, self.b, levelvec)
        points = self.grid.getPoints()
        return points

    def get_points_and_weights_component_grid(self, levelvec, numSubDiagonal):
        self.grid.setCurrentArea(self.a, self.b, levelvec)
        return self.grid.get_points_and_weights()

    def get_points_and_weights(self):
        total_points = []
        total_weights = []
        for ss in self.scheme:
            num_sub_diagonal = (self.lmax[0] + self.dim - 1) - np.sum(ss[0])
            points, weights = self.get_points_and_weights_component_grid(ss[0], num_sub_diagonal)
            total_points.extend(points)
            # adjust weights for combination -> multiply with combi coefficient
            weights = [w * ss[1] for w in weights]
            total_weights.extend(weights)
        return total_points, total_weights
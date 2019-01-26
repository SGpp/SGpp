import numpy as np
import abc, logging
from Integrator import *


# the grid class provides basic functionalities for an abstract grid
class Grid(object):

    def __init__(self, a, b, boundary=True):
        self.boundary = boundary
        self.a = a
        self.b = b

    # integrates the grid on the specified area for function f
    def integrate(self, f, levelvec, start, end):
        if not self.isGlobal():
            self.setCurrentArea(start, end, levelvec)
        return self.integrator(f, self.levelToNumPoints(levelvec), start, end)

    # def integrate_point(self, f, levelvec, start, end, point):
    #    if not self.isGlobal():
    #        self.setCurrentArea(start, end, levelvec)
    #    return self.integrator.integrate_point(f, point)

    # returns if all grid components are nested
    def isNested(self):
        return all([self.grids[d].is_nested() for d in self.dim])

    # the default case is that a grid is not globally but only locally defined
    def isGlobal(self):
        return False

    # this method can be used to generate a local grid for the given area and the specified number of points
    # typically the local grid is constructed for each area in the refinement graph during the algorithm
    @abc.abstractmethod
    def setCurrentArea(self, start, end, levelvec):
        pass

    # this method returns the actual coordinate of a point specified by the indexvector of the point
    # (e.g. in a 3 by 3 equdistant grid on the unit interval, the point [1,1] would have coordinate (0.5,0.5) and point [0,2] coordinate (0,1)
    # overwrite if necessary
    def getCoordinate(self, indexvector):
        position = np.empty(self.dim)
        for d in range(self.dim):
            position[d] = self.coordinate_array[d][indexvector[d]]
        return position

    # this method returns all the coordinate tuples of all points in the grid
    @abc.abstractmethod
    def getPoints(self):
        pass

    # this method returns the quadrature weight for the point specified by the indexvector
    @abc.abstractmethod
    def getWeight(self, indexvector):
        pass

    # this method translates a point in an equidistant mesh of level self.levelvec to its corresponding index
    def getIndexTo1DCoordinate(self, coordinate, level):
        return coordinate * 2 ** level

    def point_not_zero(self, p):
        # print(p, self.grid.boundary or not (self.point_on_boundary(p)))
        return self.boundary or not (self.point_on_boundary(p))

    def point_on_boundary(self, p):
        # print("2",p, (p == self.a).any() or (p == self.b).any())
        return (p == self.a).any() or (p == self.b).any()

    def get_points_and_weights(self):
        return self.getPoints(), self.get_weights()

    @abc.abstractmethod
    def get_weights(self):
        return list(self.getWeight(index) for index in
                    zip(*[g.ravel() for g in np.meshgrid(*[range(self.numPoints[d]) for d in range(self.dim)])]))

    def get_mid_point(self, a, b, d):
        return self.grids[d].get_mid_point(a, b)

    def setCurrentArea(self, start, end, levelvec):
        self.start = start
        self.end = end
        self.dim = len(start)
        self.levelvec = levelvec
        for d in range(self.dim):
            self.grids[d].set_current_area(self.start[d], self.end[d], self.levelvec[d])
        self.numPoints = self.levelToNumPoints(levelvec)
        self.coordinate_array = []
        self.weights = []
        self.length = np.array(end) - np.array(start)
        # prepare coordinates and weights
        for d in range(self.dim):
            self.grids[d].set_current_area(self.start[d], self.end[d], self.levelvec[d])
            coordsD, weightsD = self.grids[d].get_1d_points_and_weights()
            self.coordinate_array.append(coordsD)
            self.weights.append(weightsD)
        # print(coords)

    # this method returns the number of points in the grid that correspond to the specified levelvector
    def levelToNumPoints(self, levelvec):
        numPoints = np.zeros(len(levelvec), dtype=int)
        for d in range(len(levelvec)):
            numPoints[d] = self.grids[d].level_to_num_points_1d(levelvec[d])
        return numPoints

    def getPoints(self):
        return list(zip(*[g.ravel() for g in np.meshgrid(*self.coordinate_array)]))

    def getWeight(self, indexvector):
        weight = 1
        for d in range(self.dim):
            weight *= self.weights[d][indexvector[d]]
        return weight


from scipy.optimize import fmin
from scipy.special import eval_hermitenorm, eval_sh_legendre


class MixedGrid(Grid):
    def __init__(self, a, b, dim, grids, boundary=True, integrator=None):
        self.a = a
        self.b = b
        self.boundary = boundary
        if integrator is None:
            self.integrator = IntegratorArbitraryGridScalarProduct(self)
        else:
            if integrator == 'old':
                self.integrator = IntegratorArbitraryGrid(self)
            else:
                assert False
        self.dim = dim
        assert (len(grids) == dim)
        self.grids = grids


class Grid1d(object):
    def __init__(self, a=None, b=None, boundary=True):
        self.boundary = boundary
        self.a = a
        self.b = b

    def set_current_area(self, start, end, level):
        self.start = start
        self.end = end
        self.level = level
        self.num_points = self.level_to_num_points_1d(level)
        boundary_save = self.boundary
        self.boundary = True
        self.num_points_with_boundary = self.level_to_num_points_1d(level)
        self.boundary = boundary_save
        self.length = end - start
        # coords, weights = self.get_1d_points_and_weights(level)

        self.lowerBorder = int(0)
        self.upperBorder = self.num_points
        if not self.boundary and self.num_points < self.num_points_with_boundary:
            if start == self.a:
                self.lowerBorder = 1
            if end == self.b:
                self.upperBorder = self.num_points_with_boundary - 1
            else:
                self.upperBorder = self.num_points_with_boundary
        # equidistant spacing; only valid for equidistant grids
        self.spacing = (end - start) / (self.num_points_with_boundary - 1)

    @abc.abstractmethod
    def get_1d_points_and_weights(self):
        pass

    def get_1D_level_weights(self):
        return [self.get_1d_weight(i) for i in range(self.num_points)]

    @abc.abstractmethod
    def get_1d_weight(self, index):
        pass

    def get_mid_point(self, a, b):
        return (a + b) / 2.0

    # the default case is that a grid is nested; overwrite this if not nested!
    def is_nested(self):
        return True


# this class generates a Leja grid which constructs 1D Leja grid structures
# and constructs the tensorized grid according to the levelvector
class LejaGrid(Grid):
    def __init__(self, a, b, dim, boundary=True, integrator=None):
        self.boundary = boundary
        if integrator is None:
            self.integrator = IntegratorArbitraryGridScalarProduct(self)
        else:
            if integrator == 'old':
                self.integrator = IntegratorArbitraryGrid(self)
            else:
                assert False
        self.linear_growth_factor = 2
        self.grids = [LejaGrid1D(a=a[d], b=b[d], boundary=self.boundary) for d in range(dim)]


class LejaGrid1D(Grid1d):
    def __init__(self, a, b, boundary):
        super().__init__(a=a, b=b, boundary=boundary)
        self.linear_growth_factor = 2

    def get_1d_points_and_weights(self):
        coordsD = self.get_1D_level_points(self.level, 0, 1)
        # print(coordsD)
        weightsD = np.array(self.compute_1D_quad_weights(coordsD)) * self.length
        coordsD = np.array(coordsD)
        coordsD *= self.length
        coordsD += self.start
        return coordsD, weightsD

    def level_to_num_points_1d(self, level):
        if level == 0:
            numPoints = 1
        else:
            numPoints = self.linear_growth_factor * (level + 1) - 1
        return numPoints

    def compute_1D_quad_weights(self, grid_1D):
        N = len(grid_1D)
        V = np.zeros((N, N))

        for i in range(N):
            for j in range(N):
                V[i, j] = eval_sh_legendre(j, grid_1D[i]) * np.sqrt(2 * j + 1)

        weights = np.linalg.inv(V)[0, :]

        return weights

    def __minimize_function(self, func, a, b):

        guess = a + (b - a) / 2.0
        f_min = fmin(func, guess, xtol=1e-14, maxiter=10000, disp=False)[-1]

        return f_min

    def __get_lleja_poly(self, x, sorted_points, a, b, weightFunction=lambda x: 1.0 / s):

        if (x < a or x > b):
            return -1

        poly = 1.0
        for i in range(len(sorted_points)):
            poly *= np.abs(x - sorted_points[i])

        poly *= weightFunction(x)

        return poly

    def __get_neg_lleja_poly(self, x, sorted_points, a, b, weightFunction=lambda x: 1.0):

        return (-1.0) * self.__get_lleja_poly(x, sorted_points, a, b, weightFunction)

    def __get_starting_point(self, a, b, weightFunction=lambda x: 1.0):

        neg_weight = lambda x: -weightFunction(x)
        starting_point = self.__minimize_function(neg_weight, a, b)

        return starting_point

    def get_1D_level_points(self, curr_level, left_bound, right_bound, weightFunction=lambda x: 1.0, eps=1e-14):

        sorted_points = []
        unsorted_points = []
        no_points = self.level_to_num_points_1d(curr_level)

        starting_point = self.__get_starting_point(left_bound, right_bound, weightFunction)
        sorted_points.append(starting_point)
        unsorted_points.append(starting_point)

        for point in range(1, no_points):
            x_val = []
            y_val = []

            a = 0.
            b = left_bound
            for i in range(len(sorted_points) + 1):
                sorted_points = sorted(sorted_points)

                a = b
                if i < len(sorted_points):
                    b = sorted_points[i]
                else:
                    b = right_bound

                x_min = (a + b) / 2.0
                y_min = 0.0

                if np.abs(b - a) > eps:
                    wlleja_func = lambda x: self.__get_neg_lleja_poly(x, sorted_points, a, b, weightFunction)
                    x_min = self.__minimize_function(wlleja_func, a, b)

                    x_val.append(x_min)
                    y_min = wlleja_func(x_val[-1])
                    y_val.append(y_min)
                else:
                    x_val.append(x_min)
                    y_val.append(y_min)

            wlleja_point = x_val[y_val.index(np.min(y_val))]
            sorted_points.append(wlleja_point)
            unsorted_points.append(wlleja_point)

        unsorted_points = np.array(unsorted_points, dtype=np.float64)

        return unsorted_points


# this class provides an equdistant mesh and uses the trapezoidal rule compute the quadrature
class TrapezoidalGrid(Grid):
    def __init__(self, a, b, dim, boundary=True, integrator=None):
        self.boundary = boundary
        if integrator is None:
            self.integrator = IntegratorArbitraryGridScalarProduct(self)
        else:
            if integrator == 'old':
                self.integrator = IntegratorArbitraryGrid(self)
            else:
                assert False
        self.grids = [TrapezoidalGrid1D(a=a[d], b=b[d], boundary=self.boundary) for d in range(dim)]


class TrapezoidalGrid1D(Grid1d):

    def level_to_num_points_1d(self, level):
        return 2 ** level + 1 - (1 if not self.boundary else 0) * (
                int(1 if self.start == self.a else 0) + int(1 if self.end == self.b else 0))

    def get_1d_points_and_weights(self):
        coordsD = self.get_1D_level_points()
        weightsD = self.get_1D_level_weights()
        # print(coordsD, weightsD, self.lowerBorder, self.upperBorder, self.a, self.b, self.start, self.end)
        return coordsD, weightsD

    def get_1D_level_points(self):
        return np.linspace(self.start, self.end, self.num_points_with_boundary)[self.lowerBorder:self.upperBorder]

    def get_1d_weight(self, index):
        return self.spacing * (0.5 if index + self.lowerBorder == 0 or index + self.lowerBorder == \
                                      self.num_points_with_boundary - 1 else 1)


# this class generates a grid according to the roots of Chebyshev points and applies a Clenshaw Curtis quadrature
# the formulas are taken from: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.3141&rep=rep1&type=pdf
class ClenshawCurtisGrid(Grid):
    def __init__(self, a, b, dim, boundary=True, integrator=None):
        self.boundary = boundary
        if integrator is None:
            self.integrator = IntegratorArbitraryGridScalarProduct(self)
        else:
            if integrator == 'old':
                self.integrator = IntegratorArbitraryGrid(self)
            else:
                assert False
        self.grids = [ClenshawCurtisGrid1D(a=a[d], b=b[d], boundary=self.boundary) for d in range(dim)]


class ClenshawCurtisGrid1D(Grid1d):
    def level_to_num_points_1d(self, level):
        return 2 ** level + 1 - int(not self.boundary) * (int(self.start == 0) + int(self.end == 1))

    def get_1d_points_and_weights(self):
        coordsD = self.get_1D_level_points()
        weightsD = self.get_1D_level_weights()
        return coordsD, weightsD

    def get_1D_level_points(self):
        coordinates = np.empty(self.num_points)
        for i in range(self.num_points):
            coordinates[i] = self.start + (
                    1 - math.cos(math.pi * (i + self.lowerBorder) / (self.num_points_with_boundary - 1))) * \
                             self.length / 2
        return coordinates

    def get_1d_weight(self, index):
        weight = self.length / 2.0
        if self.num_points_with_boundary > 2:
            if index == 0 or index == self.num_points_with_boundary - 1:
                weight_factor = 1.0 / ((self.num_points_with_boundary - 2) * self.num_points_with_boundary)
            else:
                weight_factor = 0.0
                for j in range(1, math.floor((self.num_points_with_boundary - 1) / 2.0) + 1):
                    term = 1.0 / (1.0 - 4 * j * j) * math.cos(
                        2 * math.pi * index * j / (self.num_points_with_boundary - 1))
                    if j == math.floor((self.num_points_with_boundary - 1) / 2.0):
                        term *= 0.5
                    weight_factor += term
                weight_factor = 2.0 / (self.num_points_with_boundary - 1) * (1 + 2 * weight_factor)
            weight *= weight_factor
        return weight


# this class generates a grid according to the roots of Chebyshev points and applies a Clenshaw Curtis quadrature
# the formulas are taken from: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.3141&rep=rep1&type=pdf
class EquidistantGridGlobal(Grid):
    def __init__(self, a, b, boundary=True):
        self.boundary = boundary
        self.integrator = IntegratorArbitraryGrid(self)
        self.a = a
        self.b = b
        self.dim = len(a)
        self.length = np.array(b) - np.array(a)
        if (self.boundary == False):
            self.lowerBorder = np.ones(self.dim, dtype=int)
        else:
            self.lowerBorder = np.zeros(self.dim, dtype=int)

    def isGlobal(self):
        return True

    def setRefinement(self, refinement, levelvec):
        self.refinement = refinement
        self.coords = []
        self.weights = []
        self.numPoints = np.zeros(self.dim, dtype=int)
        for d in range(self.dim):
            refinementDim = refinement.get_refinement_container_for_dim(d)
            coordsD = self.mapPoints(
                refinementDim.getPointsLine()[self.lowerBorder[d]: -1 * int(self.boundary == False)], self.levelvec[d],
                d)
            weightsD = self.compute_1D_quad_weights(coords1D)
            self.coords.append(coordsD)
            self.weights.append(weightsD)
            self.numPoints[d] = len(coordsD)
        self.numPointsWithoutCoarsening = levelToNumPoints(levelvec)

    def compute_1D_quad_weights(self, grid_1D):
        N = len(grid_1D)
        V = np.zeros((N, N))

        for i in range(N):
            for j in range(N):
                V[i, j] = eval_sh_legendre(j, grid_1D[i]) * np.sqrt(2 * j + 1)

        weights = np.linalg.inv(V)[0, :]

        return weights

    def levelToNumPoints(self, levelvec):
        if hasattr(self, 'numPoints'):
            return self.numPoints
        else:
            return [2 ** levelvec[d] + 1 - int(self.boundary == False) * 2 for d in range(self.dim)]

    def getPoints(self):
        return list(zip(*[g.ravel() for g in np.meshgrid(*self.coords)]))

    def getCoordinate(self, indexvector):
        position = np.zeros(self.dim)
        for d in range(self.dim):
            position[d] = self.coords[d][indexvector[d]]
        return position

    def getWeight(self, indexvector):
        weight = 1
        for d in range(self.dim):
            weight *= self.weights[d][indexvector[d]]
        return weight

    def mapPoints(self, equidistantAdaptivePoints, level, d):
        return equidistantAdaptivePoints


# this class generates a grid according to the roots of Chebyshev points and applies a Clenshaw Curtis quadrature
# the formulas are taken from: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.3141&rep=rep1&type=pdf
class ClenshawCurtisGridGlobal(EquidistantGridGlobal):
    def __init__(self, a, b, boundary=True):
        self.boundary = boundary
        self.integrator = IntegratorArbitraryGrid(self)
        self.a = a
        self.b = b
        self.dim = len(a)
        self.length = np.array(b) - np.array(a)
        if (self.boundary == False):
            self.lowerBorder = np.ones(self.dim, dtype=int)
        else:
            self.lowerBorder = np.zeros(self.dim, dtype=int)

    '''
    def setCurrentArea(self,start,end,levelvec):
        self.start = start
        self.end = end
        self.levelvec = levelvec
        self.dim = len(self.levelvec)
        self.numPoints = self.levelToNumPoints(levelvec)
        self.numPointsWithBoundary = list(self.numPoints)
        self.length = np.array(end) - np.array(start)
        self.lowerBorder = np.zeros(self.dim,dtype=int)
        self.upperBorder = np.array(self.numPoints, dtype=int)
        if(self.boundary == False):
            for i in range(self.dim):
                if start[i]==0 :
                    self.numPointsWithBoundary[i] += 1
                    self.lowerBorder[i] = 1
                if end[i] == 1:
                    self.numPointsWithBoundary[i] += 1
                    self.upperBorder[i] = self.numPoints[i]-1
                else:
                    self.upperBorder[i] = self.numPoints[i]
    '''

    def mapPoints(self, equidistantAdaptivePoints, level, d):
        coords = np.zeros(len(equidistantAdaptivePoints))
        for i, p in enumerate(equidistantAdaptivePoints):
            index = self.getIndexTo1DCoordinate(p, level)
            coords[i] = self.get1DCoordinate(index, d)
        return coords

    def get1DCoordinate(self, index, d):
        return self.a[d] + (
                1 - math.cos(math.pi * (index + self.lowerBorder[d]) / (self.numPointsWithoutCoarsening[d] - 1))) * \
               self.length[d] / 2


import numpy.polynomial.legendre as legendre


# this class generates a grid according to the Gauss-Legendre quadrature
class GaussLegendreGrid(Grid):
    def __init__(self, a, b, dim, integrator=None):
        self.dim = dim
        self.boundary = True  # never points on boundary
        if integrator is None:
            self.integrator = IntegratorArbitraryGridScalarProduct(self)
        else:
            if integrator == 'old':
                self.integrator = IntegratorArbitraryGrid(self)
            else:
                assert False
        self.grids = [GaussLegendreGrid1D(a=a[d], b=b[d], boundary=self.boundary) for d in range(self.dim)]


class GaussLegendreGrid1D(Grid1d):
    def level_to_num_points_1d(self, level):
        return 2 ** level

    def is_nested(self):
        return False

    def get_1d_points_and_weights(self):
        coordsD, weightsD = legendre.leggauss(int(self.num_points))
        coordsD = np.array(coordsD)
        coordsD += np.ones(int(self.num_points))
        coordsD *= self.length / 2.0
        coordsD += self.start
        weightsD = np.array(weightsD) * self.length / 2
        return coordsD, weightsD


from scipy.stats import norm
from scipy.linalg import cholesky
from scipy.linalg import ldl

from scipy.sparse import diags
from scipy.sparse.linalg import lobpcg, LinearOperator
from scipy import integrate
from scipy.stats import truncnorm


# this class generates a grid according to the quadrature rule of a truncated normal distribution
# We basically compute: N * \int_a^b f(x) e^(-(x-mean)^2/(2 stddev)) dx. Where N is a normalization factor.
# The code is based on the work in "The Truncated Normal Distribution" by John Burkhardt
class TruncatedNormalDistributionGrid(Grid):
    def __init__(self, a, b, dim, mean, std_dev, integrator=None):
        # we assume here mean = 0 and std_dev = 1 for every dimension
        self.boundary = True  # never points on boundary
        self.dim = dim

        if integrator is None:
            self.integrator = IntegratorArbitraryGridScalarProduct(self)
        else:
            if integrator == 'old':
                self.integrator = IntegratorArbitraryGrid(self)
            else:
                assert False
        self.grids = [
            TruncatedNormalDistributionGrid1D(a=a[d], b=b[d], mean=mean[d], std_dev=std_dev[d], boundary=self.boundary)
            for d in range(self.dim)]
        # print(self.normalization, global_a, global_b)


class TruncatedNormalDistributionGrid1D(Grid1d):
    def __init__(self, a, b, mean, std_dev, boundary=True):
        self.shift = lambda x: x * std_dev + mean
        self.shift_back = lambda x: (x - mean) / std_dev
        self.a = self.shift_back(a)
        self.b = self.shift_back(b)
        self.normalization = 1.0 / (norm.cdf(self.b) - norm.cdf(self.a))
        self.boundary = boundary

    def get_mid_point(self, a, b):
        middle_cdf = (norm.cdf(b) + norm.cdf(a)) / 2.0
        return norm.ppf(middle_cdf)

    def level_to_num_points_1d(self, level):
        return 1 + 2 ** level

    # This grid is not nested!
    def is_nested(self):
        return False

    # this method implements the moment method for calculation the quadrature points and weights
    # a description can be found in "The Truncated Normal Distribution" from John Burkardt and
    # in "Gene Golub, John Welsch: Calculation of Gaussian Quadrature Rules"
    def get_1d_points_and_weights(self):
        num_points = int(self.num_points)
        self.L_i_dict = {}
        a = self.shift_back(self.start)
        b = self.shift_back(self.end)
        M = self.calculate_moment_matrix(num_points, a, b)

        # the stability can be improved by adding small delta to diagonal -> shifting eigenvalues away from 0
        # currently not needed
        # for i in range(num_points + 1):
        # M[i, i] += 10**-15

        # compute ldl decomposition (a bit more stable than cholesky)
        LU, diag = self.ldl_decomp(M)

        # calculate cholesky factorization based on ldl decomposition -> apply sqrt(diagonal element) to both matrices
        for i in range(num_points + 1):
            for j in range(num_points + 1):
                if (diag[j] < 0):  # in case we encounter a negative diagonal value due to instability / conditioning
                    # fill up with small value > 0
                    LU[i, j] *= 10 ** -15
                else:
                    LU[i, j] *= math.sqrt(diag[j])

        R = LU.T  # we need the upper not the lower triangular part

        # other version using scipy cholesky factorization; tends to crash earlier
        # R = scipy.linalg.cholesky(M)

        alpha_vec = np.empty(num_points)
        alpha_vec[0] = R[0, 1] / R[0, 0]
        for i in range(1, num_points):
            alpha_vec[i] = float(R[i, i + 1]) / R[i, i] - float(R[i - 1, i]) / R[i - 1, i - 1]
        beta_vec = np.empty(num_points - 1)
        for i in range(num_points - 1):
            beta_vec[i] = R[i + 1, i + 1] / R[i, i]

        # fill tridiagonal matrix

        J = np.zeros((num_points, num_points))
        for i in range(num_points):
            J[i, i] = alpha_vec[i]

        for i in range(num_points - 1):
            J[i, i + 1] = beta_vec[i]
            J[i + 1, i] = beta_vec[i]
        evals, evecs = scipy.linalg.eig(J)
        '''
        J = scipy.sparse.diags([alpha_vec, beta_vec, beta_vec], [0,-1,1])

        # get eigenvectors and eigenvalues of J
        X = np.random.rand(num_points, num_points)
        evals, evecs = scipy.sparse.linalg.lobpcg(J, X)
        '''
        points = [self.shift(ev.real) for ev in evals]
        # print(a, b, self.normalization[d])
        mu0 = self.get_moment_normalized(0, a, b)
        weights = [mu0 * value.real ** 2 for value in evecs[0]]
        # print("points and weights", num_points, a, b, points, weights)
        return points, weights

    def ldl_decomp(self, M):
        diag = np.zeros(len(M))
        L = np.zeros((len(M), len(M)))
        for i in range(len(diag)):
            L[i, i] = 1
            for j in range(i):
                summation = 0
                for k in range(j):
                    summation -= L[i, k] * L[j, k] * diag[k]
                L[i, j] = 1.0 / diag[j] * (M[i, j] + summation)
            summation = 0

            for k in range(i):
                summation -= L[i, k] ** 2 * diag[k]
            diag[i] = M[i, i] + summation
        return L, diag

    def calculate_moment_matrix(self, num_points, a, b):
        M = np.empty((num_points + 1, num_points + 1))
        for i in range(num_points + 1):
            for j in range(num_points + 1):
                M[i, j] = self.get_moment_normalized(i + j, a, b)
        return M

    # Calculation of the moments of the truncated normal distribution according to "The Truncated Normal Distribution" from John Burkardt
    # It is slightly simplified as we assume mean=0 and std_dev=1 here.
    def get_moment(self, index, a, b):
        return self.get_L_i(index, a, b)

    def get_L_i(self, index, a, b):
        if index == 0:
            return 1.0 * (norm.cdf(b) - norm.cdf(a))
        if index == 1:
            return - float((norm.pdf(b) - norm.pdf(a)))
        L_i = self.L_i_dict.get(index, None)
        if L_i is None:
            moment_m2 = self.get_L_i(index - 2, a, b)  # recursive search
            L_i = -(b ** (index - 1) * norm.pdf(b) - a ** (index - 1) * norm.pdf(a)) + (index - 1) * moment_m2
            self.L_i_dict[index] = L_i
        # print(index,L_i)
        return L_i

    # Different version of calculating moments according to: "A Recursive Formula for the Moments of a Truncated
    # Univariate Normal Distribution" by Eric Orjebin
    def get_moment2(self, index, a, b):
        if index == -1:
            return 0.0
        if index == 0:
            return 1.0 * (norm.cdf(b) - norm.cdf(a))
        moment = self.moment_dict.get(index, None)
        if moment is None:
            m_m2 = self.get_moment2(index - 2, a, b)
            moment = (index - 1) * m_m2
            moment -= b ** (index - 1) * norm.pdf(b) - a ** (index - 1) * norm.pdf(a)
            self.moment_dict[index] = moment
        return moment

    # Slight modification of 2nd version by restructuring computation; tends to be less stable
    def get_moment5(self, index, a, b):
        moment = self.moment_dict.get(index, None)
        if moment is None:
            moment = self.get_moment5_fac(index, b) * norm.pdf(b) - self.get_moment5_fac(index, a) * norm.pdf(a)
            if index % 2 == 0:
                moment += (norm.cdf(b) - norm.cdf(a))
            print(moment, index)
            self.moment_dict[index] = moment
        return moment

    def get_moment5_fac(self, index, boundary):
        if index == -1:
            return 0
        if index == 0:
            return 0
        return - (boundary ** (index - 1)) + (index - 1) * self.get_moment5_fac(index - 2, boundary)

    # Calculating moments by numerical quadrature; tends to be inaccurate
    def get_moment3(self, index, a, b):
        moment = self.moment_dict.get(index, None)
        if moment is None:
            def integrant(x):
                return x ** index * norm.pdf(x)

            moment = integrate.quad(func=integrant, a=a, b=b, epsrel=10 ** -15, epsabs=10 ** -15)[0]
            '''
            if alpha < -1 and beta > 1:
                moment = integrate.quad(func=integrant, a=alpha, b=-1, epsrel=10**-15)[0]
                moment += integrate.quad(func=integrant, a=1, b=beta, epsrel=10**-15)[0]
                moment += integrate.quad(func=integrant, a=-1, b=-0, epsrel=10**-15)[0]
                moment += integrate.quad(func=integrant, a=0, b=1, epsrel=10**-15)[0]
            '''
            # self.moment_dict[index] = moment
        return moment

    # Calculating moments using scipy
    def get_moment4(self, index, a, b):
        moment = self.moment_dict.get(index, None)
        if moment is None:
            moment = truncnorm.moment(index, a=a, b=b)
            normalization_local = 1.0 / (norm.cdf(b) - norm.cdf(a))
            moment /= normalization_local  # denormalize
            self.moment_dict[index] = moment
        return moment

    def get_moment_normalized(self, index, a, b):
        # print("Moment:", index, " variants:", self.get_moment(index,alpha,beta,mean,std_dev), self.get_moment2(index,alpha,beta,mean,std_dev), self.get_moment4(index,alpha,beta,mean,std_dev) )
        return self.get_moment(index, a, b) * self.normalization

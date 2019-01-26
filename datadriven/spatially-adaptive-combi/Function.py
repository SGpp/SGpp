import abc, logging
import numpy as np
import scipy.special
import math
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# The function class is used to define several functions for testing the algorithm
# it defines the basic interface that is used by the algorithm
class Function(object):
    # initialization if necessary
    def __init__(self):
        self.log = logging.getLogger(__name__)
        self.f_dict = {}
        self.old_f_dict = {}
        self.do_cache = True  # indicates whether function values should be cached

    def reset_dictionary(self):
        self.old_f_dict = {**self.old_f_dict, **self.f_dict}
        self.f_dict = {}

    def __call__(self, coordinates):
        if self.do_cache:
            coords = tuple(coordinates)
            f_value = self.f_dict.get(coords, None)
            if f_value is not None:
                return f_value
            f_value = self.old_f_dict.get(coords, None)
            if f_value is not None:
                self.f_dict[coords] = f_value
                return f_value
        f_value = self.eval(coords)
        if self.do_cache:
            self.f_dict[coords] = f_value
        return f_value

    def deactivate_caching(self):
        self.do_cache = False

    def get_f_dict_size(self):
        return len(self.f_dict)

    def get_f_dict_points(self):
        return list(self.f_dict.keys())

    def get_f_dict_values(self):
        return list(self.f_dict.values())

    # evaluates the function at the specified coordinate
    @abc.abstractmethod
    def eval(self, coordinates):
        return

    # this returns the analytic solution of the integral in the specified area
    # currently necessary for the error estimator
    @abc.abstractmethod
    def getAnalyticSolutionIntegral(self, start, end):
        return

    # this method plots the function in the specified area for 2D
    def plot(self, start, end, filename=None):
        dim = len(start)
        if dim > 2:
            print("Cannot plot function with dim > 2")
            return
        xArray = np.linspace(start[0], end[0], 10 ** 2)
        yArray = np.linspace(start[1], end[1], 10 ** 2)
        X = [x for x in xArray]
        Y = [y for y in yArray]
        X, Y = np.meshgrid(X, Y)
        Z = np.zeros(np.shape(X))
        for i in range(len(X)):
            for j in range(len(X[i])):
                # print(X[i,j],Y[i,j],self.eval((X[i,j],Y[i,j])))
                Z[i, j] = self.eval((X[i, j], Y[i, j]))
        # Z=self.eval((X,Y))
        # print Z
        fig = plt.figure(figsize=(14, 6))

        # `ax` is a 3D-aware axis instance, because of the projection='3d' keyword argument to add_subplot
        ax = fig.add_subplot(1, 2, 1, projection='3d')

        # p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False)
        # plt.show()
        if filename is not None:
            fig.savefig(filename, bbox_inches='tight')
        plt.show()


from scipy import integrate
from scipy.stats import norm


class FunctionShift(Function):
    def __init__(self, function, shift):
        super().__init__()
        self.function = function
        self.shift = shift  # a list of functions that can manipulate the coordinates for every dimension

    def eval(self, coordinates):
        shifted_coordinates = self.shift(coordinates)
        return self.function.eval(shifted_coordinates)

    def getAnalyticSolutionIntegral(self, start, end):
        start_shifted = self.shift(start)
        end_shifted = self.shift(end)
        return self.function.getAnalyticSolutionIntegral(start_shifted, end_shifted)


class FunctionUQNormal(Function):
    def __init__(self, function, mean, std_dev, a, b):
        super().__init__()
        self.dim = len(mean)
        self.mean = mean
        self.std_dev = std_dev
        self.function = function
        self.a_global = a
        self.b_global = b

    def eval(self, coordinates):
        return self.function(np.asarray(coordinates) * np.asarray(self.std_dev) + np.asarray(self.mean))

    def eval_with_normal(self, coordinates):
        value = self.eval(coordinates)
        dim = len(coordinates)
        # add contribution of normal distribution
        summation = 0
        for d in range(dim):
            summation -= (coordinates[d]) ** 2 / (2 * 1)
        return value * np.exp(summation)

    def getAnalyticSolutionIntegral(self, start, end):
        if self.dim == 3:
            f = lambda x, y, z: self.eval_with_normal([x, y, z])
            normalization = 1
            for d in range(len(start)):
                S = norm.cdf(self.b_global[d]) - norm.cdf(self.a_global[d])
                normalization *= 1.0 / (S * math.sqrt(2 * math.pi * 1))
            return normalization * \
                   integrate.tplquad(f, start[2], end[2], lambda x: start[1], lambda x: end[1], lambda x, y: start[0],
                                     lambda x, y: end[0])[0]
        elif self.dim == 2:
            f = lambda x, y: self.eval_with_normal([x, y])
            normalization = 1
            for d in range(len(start)):
                S = norm.cdf(self.b_global[d]) - norm.cdf(self.a_global[d])
                normalization *= 1.0 / (S * math.sqrt(2 * math.pi * 1))
            return normalization * integrate.dblquad(f, start[1], end[1], lambda x: start[0], lambda x: end[0])[0]


class FunctionUQNormal2(Function):
    def __init__(self, function, mean, std_dev, a, b):
        super().__init__()
        self.mean = mean
        self.std_dev = std_dev
        self.function = function
        self.a_global = a
        self.b_global = b
        self.dim = len(mean)

    def eval(self, coordinates):
        return self.function.eval(coordinates)

    def eval_with_normal(self, coordinates):
        dim = len(coordinates)
        # add contribution of normal distribution
        return np.prod([norm.pdf(x=coordinates[d], loc=self.mean[d], scale=self.std_dev[d]) / (
                    norm.cdf(self.b_global[d], loc=self.mean[d], scale=self.std_dev[d]) - norm.cdf(self.a_global[d],
                                                                                                   loc=self.mean[d],
                                                                                                   scale=self.std_dev[
                                                                                                       d])) for d in
                        range(dim)])

    def getAnalyticSolutionIntegral(self, start, end):
        if self.dim == 3:
            f = lambda x, y, z: self.eval([x, y, z]) * self.eval_with_normal([x, y, z])
            return integrate.tplquad(f, start[2], end[2], lambda x: start[1], lambda x: end[1], lambda x, y: start[0],
                                     lambda x, y: end[0])[0]
        elif self.dim == 2:
            f = lambda x, y: self.eval([x, y]) * self.eval_with_normal([x, y])
            return integrate.dblquad(f, start[1], end[1], lambda x: start[0], lambda x: end[0])[0]
        else:
            assert False


class FunctionUQWeighted(Function):
    def __init__(self, function, weight_function):
        super().__init__()
        self.function = function
        self.weight_function = weight_function

    def eval(self, coordinates):
        return self.function.eval(coordinates)

    def getAnalyticSolutionIntegral(self, start, end):
        self.dim = len(start)
        if self.dim == 3:
            f = lambda x, y, z: self.eval([x, y, z]) * self.weight_function([x, y, z])
            return \
                integrate.tplquad(f, start[2], end[2], lambda x: start[1], lambda x: end[1], lambda x, y: start[0],
                                  lambda x, y: end[0])[0]
        elif self.dim == 2:
            f = lambda x, y: self.eval([x, y]) * self.weight_function([x, y])
            return integrate.dblquad(f, start[1], end[1], lambda x: start[0], lambda x: end[0])[0]
        else:
            assert False


class FunctionUQ(Function):
    def eval(self, coordinates):
        assert (len(coordinates) == 3)
        parameter1 = coordinates[0]
        parameter2 = coordinates[1]
        parameter3 = coordinates[2]

        # Model with discontinuity
        # Nicholas Zabarras Paper: „Sparse grid collocation schemes for stochastic natural convection problems“
        # e^(-x^2 + 2*sign(y))
        value_of_interest = math.exp(-parameter1 ** 2 + 2 * np.sign(parameter2)) + parameter3
        return value_of_interest

    def getAnalyticSolutionIntegral(self, start, end):
        f = lambda x, y, z: self.eval([x, y, z])
        return integrate.tplquad(f, start[2], end[2], lambda x: start[1], lambda x: end[1], lambda x, y: start[0],
                                 lambda x, y: end[0])[0]


from scipy.stats import truncnorm


class FunctionUQ2(Function):
    def eval(self, coordinates):
        # print(coordinates)
        assert (len(coordinates) == 2)
        parameter1 = coordinates[0]
        parameter2 = coordinates[1]

        # Model with discontinuity
        # Nicholas Zabarras Paper: „Sparse grid collocation schemes for stochastic natural convection problems“
        # e^(-x^2 + 2*sign(y))
        value_of_interest = math.exp(-parameter1 ** 2 + 2 * np.sign(parameter2))
        return value_of_interest

    def getAnalyticSolutionIntegral(self, start, end):
        f = lambda x, y: self.eval([x, y])
        return integrate.dblquad(f, start[1], end[1], lambda x: start[0],
                                 lambda x: end[0])[0]


# This class composes different functions that fullfill the function interface to generate a composition of functions
# Each Function can be weighted by a factor to allow for a flexible composition
class FunctionCompose(Function):
    def __init__(self, functions):
        super().__init__()
        self.functions = functions

    def eval(self, coordinates):
        result = 0.0
        for (f, factor) in self.functions:
            result += f.eval(coordinates) * factor
        return result

    def getAnalyticSolutionIntegral(self, start, end):
        result = 0.0
        for (f, factor) in self.functions:
            result += f.getAnalyticSolutionIntegral(start, end) * factor
        return result


# This Function represents the corner Peak f the genz test functions
class GenzCornerPeak(Function):
    def __init__(self, coeffs):
        super().__init__()
        self.coeffs = np.array(coeffs)
        self.dim = len(coeffs)

    def eval(self, coordinates):
        result = 1
        for d in range(self.dim):
            result += self.coeffs[d] * coordinates[d]
        return result ** (-self.dim - 1)

    def getAnalyticSolutionIntegral(self, start, end):
        factor = ((-1) ** self.dim) * 1.0 / (math.factorial(self.dim) * np.prod(self.coeffs))
        combinations = list(zip(*[g.ravel() for g in np.meshgrid(*[[0, 1] for d in range(self.dim)])]))
        result = 0
        for c in combinations:
            partial_result = 1
            for d in range(self.dim):
                if c[d] == 1:
                    value = start[d]
                else:
                    value = end[d]
                partial_result += value * self.coeffs[d]
            result += (-1) ** sum(c) * partial_result ** -1
        return factor * result


class GenzProductPeak(Function):
    def __init__(self, coefficients, midpoint):
        super().__init__()
        self.coeffs = coefficients
        self.midPoint = midpoint
        self.dim = len(coefficients)
        self.factor = 10 ** (-self.dim)

    def eval(self, coordinates):
        result = 1
        for d in range(self.dim):
            result /= (self.coeffs[d] ** (-2) + (coordinates[d] - self.midPoint[d]) ** 2)
        return result * self.factor

    def getAnalyticSolutionIntegral(self, start, end):
        result = 1
        for d in range(self.dim):
            result *= np.arctan(self.coeffs[d] * (self.midPoint[d] - start[d])) * self.coeffs[d] - np.arctan(
                self.coeffs[d] * (self.midPoint[d] - end[d])) * self.coeffs[d]
        return result * self.factor


class GenzOszillatory(Function):
    def __init__(self, coeffs, offset):
        super().__init__()
        self.coeffs = coeffs
        self.dim = len(coeffs)
        self.offset = offset

    def eval(self, coordinates):
        result = 2 * math.pi * self.offset
        for d in range(self.dim):
            result += self.coeffs[d] * coordinates[d]
        return math.cos(result)

    def getAnalyticSolutionIntegral(self, start, end):
        factor = ((-1) ** int(math.floor(self.dim / 2))) * 1.0 / np.prod(self.coeffs)
        combinations = list(zip(*[g.ravel() for g in np.meshgrid(*[[0, 1] for d in range(self.dim)])]))
        result = 0
        for c in combinations:
            partial_result = 2 * math.pi * self.offset
            for d in range(self.dim):
                if c[d] == 1:
                    value = start[d]
                else:
                    value = end[d]
                partial_result += value * self.coeffs[d]
            if self.dim % 2 == 1:
                result += (-1) ** sum(c) * math.sin(partial_result)
            else:
                result += (-1) ** sum(c) * math.cos(partial_result)
        return factor * result


class GenzDiscontinious(Function):
    def __init__(self, coeffs, border):
        super().__init__()
        self.coeffs = coeffs
        self.border = border
        self.dim = len(coeffs)

    def eval(self, coordinates):
        result = 0
        for d in range(self.dim):
            if coordinates[d] >= self.border[d]:
                return 0.0
            result -= self.coeffs[d] * coordinates[d]
        return np.exp(result)

    def getAnalyticSolutionIntegral(self, start, end):
        result = 1
        end = list(end)
        for d in range(self.dim):
            if start[d] >= self.border[d]:
                return 0.0
            else:
                end[d] = min(end[d], self.border[d])
                result *= (np.exp(-self.coeffs[d] * start[d]) - np.exp(-self.coeffs[d] * end[d])) / self.coeffs[d]
        return result


class GenzC0(Function):
    def __init__(self, coeffs, midpoint):
        super().__init__()
        self.coeffs = coeffs
        self.midPoint = midpoint
        self.dim = len(coeffs)

    def eval(self, coordinates):
        result = 0
        for d in range(self.dim):
            result -= self.coeffs[d] * abs(coordinates[d] - self.midPoint[d])
        return np.exp(result)

    def getAnalyticSolutionIntegral(self, start, end):
        result = 1
        for d in range(self.dim):
            one_d_integral = 0
            if start[d] < self.midPoint[d]:
                if end[d] < self.midPoint[d]:
                    one_d_integral += (np.exp(self.coeffs[d] * (end[d] - self.midPoint[d]))) / self.coeffs[d] - (
                        np.exp(self.coeffs[d] * (start[d] - self.midPoint[d]))) / self.coeffs[d]
                else:
                    one_d_integral += 1.0 / self.coeffs[d] - (
                        np.exp(self.coeffs[d] * (start[d] - self.midPoint[d]))) / self.coeffs[d]
            if end[d] > self.midPoint[d]:
                if start[d] > self.midPoint[d]:
                    one_d_integral += (np.exp(self.coeffs[d] * (self.midPoint[d] - start[d]))) / self.coeffs[d] - (
                        np.exp(self.coeffs[d] * (self.midPoint[d] - end[d]))) / self.coeffs[d]
                else:
                    one_d_integral += 1.0 / self.coeffs[d] - (
                        np.exp(self.coeffs[d] * (self.midPoint[d] - end[d]))) / self.coeffs[d]
            result *= one_d_integral
        return result


# This function is the test case function 2 of the paper from Jakeman and Roberts: https://arxiv.org/pdf/1110.0010.pdf
class Function2Jakeman(Function):
    def __init__(self, midpoint, coefficients):
        super().__init__()
        self.midpoint = midpoint
        self.coefficients = coefficients

    def eval(self, coordinates):
        dim = len(coordinates)
        assert (dim == len(self.coefficients))
        summation = 0.0
        for d in range(dim):
            summation -= self.coefficients[d] * (coordinates[d] - self.midpoint[d]) ** 2
        return np.exp(summation)

    def getAnalyticSolutionIntegral(self, start, end):
        dim = len(start)
        # print lowerBounds,upperBounds,coefficients, midpoints
        result = 1.0
        sqPiHalve = np.sqrt(np.pi) * 0.5
        for d in range(dim):
            result = result * (
                    sqPiHalve * scipy.special.erf(np.sqrt(self.coefficients[d]) * (end[d] - self.midpoint[d])) -
                    sqPiHalve * scipy.special.erf(
                np.sqrt(self.coefficients[d]) * (start[d] - self.midpoint[d]))) / np.sqrt(self.coefficients[d])
        return result


# This function is the first test function of the paper from Gerstner and Griebel:
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.3141&rep=rep1&type=pdf
class FunctionGriebel(Function):
    def eval(self, coordinates):
        dim = len(coordinates)
        prod = 1.0
        for d in range(dim):
            prod *= coordinates[d] ** (1.0 / dim)
        return (1 + 1.0 / dim) ** dim * prod

    def getAnalyticSolutionIntegral(self, start, end):
        dim = len(start)
        result = 1.0
        for d in range(dim):
            result = result * (
                    end[d] ** (1 + 1.0 / dim) / (1 + 1.0 / dim) - start[d] ** (1 + 1.0 / dim) / (1 + 1.0 / dim))
        return (1 + 1.0 / dim) ** dim * result


# This is a variant of the generalized normal distribution
# Currently the analytic solution is not correct!!!
class FunctionGeneralizedNormal(Function):
    def __init__(self, midpoints, coefficients, exp):
        super().__init__()
        self.midpoints = midpoints
        self.coefficients = coefficients
        self.exponent = exp

    def eval(self, coordinates):
        dim = len(coordinates)
        assert (dim == len(self.coefficients))
        summation = 0.0
        for d in range(dim):
            summation -= self.coefficients[d] * (abs(coordinates[d] - self.midpoints[d])) ** 2
        return np.exp(summation ** self.exponent)

    def getAnalyticSolutionIntegral(self, start, end):
        dim = len(start)
        # print lowerBounds,upperBounds,coefficients, midpoints
        result = 1.0
        sqPiHalve = np.sqrt(np.pi) * 0.5
        for d in range(dim):
            result = result * (
                    sqPiHalve * scipy.special.erf(np.sqrt(self.coefficients[d]) * (end[d] - self.midpoints[d])) -
                    sqPiHalve * scipy.special.erf(
                np.sqrt(self.coefficients[d]) * (start[d] - self.midpoints[d]))) / np.sqrt(self.coefficients[d])
        return result

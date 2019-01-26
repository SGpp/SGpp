import scipy.integrate
import numpy as np
import abc
import logging


# This class is the general interface of an error estimator currently used by the algorithm
class ErrorCalculator(object):
    # initialization
    def __init__(self):
        self.log = logging.getLogger(__name__)

    # calculates error for the function f and the integral information that was computed by the algorithm
    # this information contians the area specification and the approximated integral
    # current form is (approxIntegral,start,end)
    @abc.abstractmethod
    def calc_error(self, f, refine_object):
        return


# This error estimator doea a comparison to analytic solution. It outputs the absolute error.
class ErrorCalculatorAnalytic(ErrorCalculator):
    def calc_error(self, f, refine_object):
        lower_bounds = refine_object.start
        upper_bounds = refine_object.end
        real_integral_value = f.getAnalyticSolutionIntegral(lower_bounds, upper_bounds)
        return abs(refine_object.integral - real_integral_value)


# This error estimator doea a comparison to analytic solution. It outputs the relative error.
class ErrorCalculatorAnalyticRelative(ErrorCalculator):
    def calc_error(self, f, refine_object):
        lower_bounds = refine_object.start
        upper_bounds = refine_object.end
        real_integral_value = f.getAnalyticSolutionIntegral(lower_bounds, upper_bounds)
        real_integral_complete = f.getAnalyticSolutionIntegral(a, b)
        return abs((refine_object.integral - real_integral_value) / real_integral_complete)


# This error estimator does a surplus estimation. It outputs the absolute error.
class ErrorCalculatorSurplusCell(ErrorCalculator):
    def calc_error(self, f, refine_object):
        error = self.calc_area_error(refine_object.sub_integrals)
        return error

    def calc_area_error(self, sub_integrals):
        error = 0.0
        for sub_integral in sub_integrals:
            error += sub_integral[0] * sub_integral[1]
        return abs(error)


class ErrorCalculatorSurplusCellPunishDepth(ErrorCalculatorSurplusCell):
    def calc_error(self, f, refine_object):
        lower_bounds = np.array(refine_object.start)
        upper_bounds = np.array(refine_object.end)
        error = self.calc_area_error(refine_object.sub_integrals)
        return error * np.prod(upper_bounds - lower_bounds)
        

class ErrorCalculatorExtendSplit(ErrorCalculator): #TODO
    def calc_error(self, f, refine_object):
        return abs(refine_object.integral - refine_object.parent_integral)


class ErrorCalculatorExtendSplit2(ErrorCalculator):  # TODO
    def calc_error(self, f, refine_object):
        return abs(refine_object.integral - refine_object.parent_integral) / refine_object.factor
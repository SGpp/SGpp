import abc, logging
# Python modules
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import operator
import numpy as np
import scipy as sp
import scipy.integrate
from scipy.interpolate import interpn
import scipy.special
import math
import time
from RefinementContainer import *
from RefinementObject import *
from combiScheme import *
from Grid import *
from ErrorCalculator import *
from Function import *
from StandardCombi import *


# This class defines the general interface and functionalties of all spatially adaptive refinement strategies
class SpatiallyAdaptivBase(StandardCombi):
    def __init__(self, a, b, grid=None):
        self.log = logging.getLogger(__name__)
        self.dim = len(a)
        self.a = a
        self.b = b
        self.grid = grid
        self.refinements_for_recalculate = 100
        assert (len(a) == len(b))

    # returns the number of points in a single component grid with refinement
    def get_num_points_component_grid(self, levelvec, do_naive, num_sub_diagonal):
        array2 = self.get_points_component_grid(levelvec, num_sub_diagonal)
        if do_naive:
            array2new = array2
        else:  # remove points that appear in the list multiple times
            array2new = list(set(array2))
        # print(len(array2new))
        return len(array2new)

    def evaluate_final_combi(self):
        combiintegral = 0
        dim = self.dim
        # print "Dim:",dim
        num_evaluations = 0
        for ss in self.scheme:
            integral = 0
            for area in self.get_areas():
                area_integral, partial_integrals, evaluations = self.evaluate_area(self.f, area, ss[0])
                if area_integral != -2 ** 30:
                    num_evaluations += evaluations
                    integral += area_integral
            integral *= ss[1]
            combiintegral += integral
        return combiintegral, num_evaluations

    def init_adaptive_combi(self, f, minv, maxv, refinement_container, tol):
        self.tolerance = tol
        self.f = f
        if self.realIntegral is not None:
            print("Reference solution:", self.realIntegral)
        else:
            print("No reference solution present. Working purely on surplus error estimates.")
        if (refinement_container == []):  # initialize refinement
            self.lmin = [minv for i in range(self.dim)]
            self.lmax = [maxv for i in range(self.dim)]
            # calculate the combination scheme
            self.combischeme = CombiScheme(self.dim)
            self.scheme = self.combischeme.getCombiScheme(self.lmin[0], self.lmax[0], self.dim)
            self.initialize_refinement()
            self.f.reset_dictionary()
        else:  # use the given refinement; in this case reuse old lmin and lmax and finestWidth; works only if there was no other run in between on same object
            self.refinement = refinement_container
            self.refinement.reinit_new_objects()
        # initialize values
        self.refinements = 0
        # self.combiintegral = 0
        # self.subAreaIntegrals = []
        self.counter = 1
        # self.evaluationsTotal = 0 #number of evaluations in current grid
        # self.evaluationPerArea = [] #number of evaluations per area

    def evaluate_integral(self):
        # initialize values
        # number_of_evaluations = 0
        # get tuples of all the combinations of refinement to access each subarea (this is the same for each component grid)
        areas = self.get_new_areas()
        integralarrayComplete = np.zeros(len(areas))
        evaluation_array = np.zeros(len(areas))
        # calculate integrals
        for ss in self.scheme:  # iterate over component grids
            # initialize component grid specific variables
            numSubDiagonal = (self.lmax[0] + self.dim - 1) - np.sum(ss[0])
            integral = 0
            # iterate over all areas and calculate the integral
            for k, area in enumerate(areas):
                # print(ss)
                area_integral, partial_integrals, evaluations = self.evaluate_area(self.f, area, ss[0])
                if area_integral != -2 ** 30:
                    if partial_integrals is not None:  # outdated
                        pass
                        # integralArrayIndividual.extend(partial_integrals)
                    else:
                        integralarrayComplete[k] += ss[1] * area_integral
                        # self.combiintegral += area_integral * ss[1]
                        evaluation_array[k] += evaluations

        for k in range(len(integralarrayComplete)):
            i = k + self.refinement.size() - self.refinement.new_objects_size()
            self.refinement.set_integral(i, integralarrayComplete[k])
            self.refinement.set_evaluations(i, evaluation_array[k] / len(self.scheme))
            self.calc_error(i, self.f)

        # getArea with maximal error
        self.errorMax = self.refinement.get_max_error()
        self.total_error = self.refinement.get_total_error()
        print("max surplus error:", self.errorMax, "total surplus error:", self.total_error)
        if self.realIntegral is not None:
            return abs(self.refinement.integral - self.realIntegral) / abs(self.realIntegral)
        else:
            return self.total_error

    def refine(self):
        # split all cells that have an error close to the max error
        areas = self.get_areas()
        self.prepare_refinement()
        self.refinement.clear_new_objects()
        margin = 0.9
        quit_refinement = False
        while True:  # refine all areas for which area is within margin
            # get next area that should be refined
            found_object, position, refine_object = self.refinement.get_next_object_for_refinement(
                tolerance=self.errorMax * margin)
            if found_object and not quit_refinement:  # new area found for refinement
                self.refinements += 1
                # print("Refining position", position)
                quit_refinement = self.do_refinement(refine_object, position)

            else:  # all refinements done for this iteration -> reevaluate integral and check if further refinements necessary
                print("Finished refinement")
                self.refinement_postprocessing()
                break

        if self.recalculate_frequently and self.refinements / self.refinements_for_recalculate > self.counter:
            self.refinement.reinit_new_objects()
            self.combiintegral = 0
            self.subAreaIntegrals = []
            self.evaluationPerArea = []
            self.evaluationsTotal = 0
            self.counter += 1
            print("recalculating errors")

    # optimized adaptive refinement refine multiple cells in close range around max variance (here set to 10%)
    def performSpatiallyAdaptiv(self, minv=1, maxv=2, f=FunctionGriebel(), errorOperator=None, tol=10 ** -2,
                                refinement_container=[], do_plot=False, recalculate_frequently=False, test_scheme=False,
                                reevaluate_at_end=False, reference_solution=None):
        self.errorEstimator = errorOperator
        self.recalculate_frequently = recalculate_frequently
        self.realIntegral = reference_solution

        self.init_adaptive_combi(f, minv, maxv, refinement_container, tol)
        error_array = []
        num_point_array = []
        while True:
            error = self.evaluate_integral()
            error_array.append(error)
            num_point_array.append(self.get_total_num_points(distinct_function_evals=True))
            print("combiintegral:", self.refinement.integral)
            print("Current error:", error)
            # check if tolerance is already fullfilled with current refinement
            if error > tol:
                # refine further
                self.refine()
                if do_plot:
                    print("Refinement Graph:")
                    self.draw_refinement()
                    print("Combi Scheme:")
                    self.print_resulting_combi_scheme()
                    print("Resulting Sparse Grid:")
                    self.print_resulting_sparsegrid()
            else:  # refinement finished
                break
        # finished adaptive algorithm
        print("Number of refinements", self.refinements)
        if test_scheme:
            self.check_combi_scheme()
        if reevaluate_at_end:
            # evaluate final integral
            combiintegral, number_of_evaluations = self.evaluate_final_combi()
        else:
            combiintegral = self.refinement.integral
            number_of_evaluations = self.refinement.evaluationstotal
        return self.refinement, self.scheme, self.lmax, combiintegral, number_of_evaluations, error_array, num_point_array

    @abc.abstractmethod
    def initialize_refinement(self):
        pass

    @abc.abstractmethod
    def get_points_component_grid(self, levelvec, numSubDiagonal):
        return

    @abc.abstractmethod
    def evaluate_area(self, f, area, levelvec):
        pass

    @abc.abstractmethod
    def do_refinement(self, area, position):
        pass

    # this is a default implementation that should be overritten if necessary
    def prepare_refinement(self):
        pass

    # this is a default implementation that should be overritten if necessary
    def refinement_postprocessing(self):
        self.refinement.apply_remove()
        self.refinement.refinement_postprocessing()

    # this is a default implementation that should be overritten if necessary
    def calc_error(self, objectID, f):
        self.refinement.calc_error(objectID, f)

    # this is a default implementation that should be overritten if necessary
    def get_new_areas(self):
        return self.refinement.get_new_objects()

    # this is a default implementation that should be overritten if necessary
    def get_areas(self):
        return self.refinement.get_objects()

    # this method can be overwritten if for the method a graphical refinement visualization exists
    def draw_refinement(self, filename=None):
        pass

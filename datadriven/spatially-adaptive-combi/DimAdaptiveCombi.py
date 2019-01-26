from StandardCombi import *
from combiScheme import *
from Grid import *


# T his class implements the standard combination technique
class DimAdaptiveCombi(StandardCombi):
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

    # standard combination scheme for quadrature
    # lmin = minimum level; lmax = target level
    # f = function to integrate; dim=dimension of problem
    def perform_combi(self, minv, maxv, f, tolerance, reference_solution=None):
        start = self.a
        end = self.b
        self.f = f
        self.f.reset_dictionary()
        # compute minimum and target level vector
        self.lmin = [minv for i in range(self.dim)]
        self.lmax = [maxv for i in range(self.dim)]
        real_integral = reference_solution
        assert(reference_solution is not None)
        self.combischeme.init_adaptive_combi_scheme(maxv, minv)
        combiintegral = 0
        self.scheme = self.combischeme.getCombiScheme(self.lmin[0], self.lmax[0], self.dim)
        integral_dict = {}
        errors = []  # tracks the error evolution during the refinement procedure
        num_points = []  # tracks the number of points during the refinement procedure
        while True:
            combiintegral = 0
            self.scheme = self.combischeme.getCombiScheme(self.lmin[0], self.lmax[0], self.dim, do_print=False)
            error_array = np.zeros(len(self.scheme))
            for i, ss in enumerate(self.scheme):
                if tuple(ss[0]) not in integral_dict:
                    integral = self.grid.integrate(self.f, ss[0], start, end)
                    integral_dict[tuple(ss[0])] = integral
                else:
                    integral = integral_dict[tuple(tuple(ss[0]))]
                # as error estimator we compare to the analytic solution and divide by the cost=number of points in grid
                error_array[i] = abs(integral - real_integral) / abs(real_integral) / np.prod(
                    self.grid.levelToNumPoints(ss[0])) if self.combischeme.is_refinable(ss[0]) else 0
                combiintegral += integral * ss[1]
            do_refine = True
            if abs(combiintegral - real_integral) / abs(real_integral) < tolerance:
                break
            print("Current combi integral:", combiintegral)
            print("Currentrelative error:", abs(combiintegral - real_integral) / abs(real_integral))
            errors.append(abs(combiintegral - real_integral) / abs(real_integral))
            num_points.append(self.get_total_num_points(distinct_function_evals=True))
            while do_refine:
                grid_id = np.argmax(error_array)
                # print(error_array)
                print("Current error:", abs(combiintegral - real_integral) / abs(real_integral))
                print("Refining", self.scheme[grid_id], self.combischeme.is_refinable(self.scheme[grid_id][0]))
                refined_dims = self.combischeme.update_adaptive_combi(self.scheme[grid_id][0])
                do_refine = refined_dims == []
                error_array[grid_id] = 0.0
        print("Final scheme:")
        self.scheme = self.combischeme.getCombiScheme(self.lmin[0], self.lmax[0], self.dim, do_print=True)
        print("CombiSolution", combiintegral)
        print("Analytic Solution", real_integral)
        print("Difference", abs(combiintegral - real_integral))
        return self.scheme, abs(combiintegral - real_integral), combiintegral, errors, num_points

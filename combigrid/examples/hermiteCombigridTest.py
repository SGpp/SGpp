import pysgpp
import FunctionClasses as fctClass
import matplotlib.pyplot as plt
import plotlib as p
import numpy as np
import GridCollection as gC
import random
import saveLoad as sl

import sympy as sp
from mpl_toolkits.mplot3d import Axes3D


def operationwrapper(combi_operation, level):
    def operation(x):
        return combi_operation.evaluate(level, x)

    return operation


# for functions from sgpp::optimization
def getfuncwrapper(func):
    def function(x):
        return func.evalUndisplaced(x)

    return function


def estimatel2Error(n, dim, gridOpEval, targetFunc):
    sum = 0

    for _ in range(n):
        point = pysgpp.DataVector(dim)
        for d in range(dim):
            point[d] = random.random()
        #print(gridOpEval(point))
        sum += (gridOpEval(point) - targetFunc(point)) ** 2
    sum = sum / n
    sum = np.sqrt(sum)
    return sum


def estimatel2ErrorGradients(n, dim, gridOpEval_container, func_container):
    """

    Returns:
        List: the l2 error in the gradients as a list in powerset order
    """

    mixed_directions = list(fctClass.powerset([i for i in range(dim)]))

    mixed_grad = []
    opEval_mixgrad = []

    # do the derivatives for operation, function and save resulting functions in the 2 arrays
    for i in range(1, len(mixed_directions)):
        opEval_mixgrad_temp = gridOpEval_container.getGradient(mixed_directions[i])
        mixed_grad_temp = func_container.getGradient(mixed_directions[i])
        mixed_grad.append(pysgpp.multiFunc(mixed_grad_temp))
        opEval_mixgrad.append(opEval_mixgrad_temp)

    errors = []
    for i in range(len(mixed_grad)):
        errors.append(estimatel2Error(n, dim, opEval_mixgrad[i], mixed_grad[i]))

    return errors


def calc_error_on_gridpts(gridpoints, gridOpEval, targetFunc):
    error = 0

    for i in range(len(gridpoints)):
        x = gridpoints[i]
        error = max(gridOpEval(x) - targetFunc(x), error)
        # print(str(x[0]) + "   " + str(x[1]), abs(gridOpEval(x) - targetFunc(x)), gridOpEval(x),
        # targetFunc(x))

    return error


def calc_error_gradient(gridpoints, gridOpEval, targetFunc, dim):
    error = []

    # iterate over dimensions
    for d in range(dim):
        grad_i = fctClass.getgradkfunc(targetFunc, d)
        opEvalgrad_i = fctClass.getgradkfunc(gridOpEval, d)
        error.append(calc_error_on_gridpts(gridpoints, opEvalgrad_i, grad_i))
        # print("-------------------------------------------------------------")

    return error


# calculates the error from gradients on the gridpoints
def calc_error_mixed_gradient(gridpoints, gridOpEval_container, func_container, dim):
    # type: (object, fctClass.funcGradientCollection, fctClass.funcGradientCollection, int) ->list

    mixed_directions = list(fctClass.powerset([i for i in range(dim)]))

    mixed_grad = []
    opEval_mixgrad = []

    for i in range(1, len(mixed_directions)):
        opEval_mixgrad_temp = gridOpEval_container.getGradient(mixed_directions[i])
        func_mixed_grad_temp = func_container.getGradient(mixed_directions[i])
        mixed_grad.append(pysgpp.multiFunc(func_mixed_grad_temp))
        opEval_mixgrad.append(opEval_mixgrad_temp)

    errors = []
    for i in range(len(mixed_grad)):
        errors.append(calc_error_on_gridpts(gridpoints, opEval_mixgrad[i], mixed_grad[i]))

    return errors


def calc_tangent(x, k, operation):
    grad_k = fctClass.getgradkfunc(operation, k)

    m = grad_k(x)
    b = operation(x)
    return m, b


# returns a number of points to draw the tangent
def tangent_samples(x, dim, k, m, b):
    samples = 2
    x_samples = []
    for d in range(dim):
        x_samples.append(np.full(samples, x[d]))

    x_samples[k] = (np.linspace(x[k] - 0.05, x[k] + 0.05, samples))

    y = [(x_samples[k][i] - x[k]) * m + b for i in range(samples)]

    return x_samples[0], x_samples[1], y


def calc_gradient_tangent_gridpoint(gridpoints, operation, dim):
    x0_all = []
    x1_all = []
    y_all = []
    for p in range(len(gridpoints)):

        for d in range(dim):
            m, b = calc_tangent(gridpoints[p], d, operation)

            x0, x1, y = tangent_samples(gridpoints[p], 2, d, m, b)
            x0_all.append(x0)
            x1_all.append(x1)
            y_all.append(y)

    return x0_all, x1_all, y_all


# parabola between [0,1]
def f1D(x):
    #return -4 * ((x[0] - 0.5) ** 2) + 1
    return x[0]


def f1D_grad(x):
    return 4 - 8 * x[0]


def f2D(x):
    return 4 * x[0] * (1 - x[0]) * 4 * x[1] * (1 - x[1])


def f2D_test(x):
    return x[0] ** 2 + x[1] ** 2


# returns the plotgridpoints and the points wrapped as DataVector
def generate1DGrid(samples):
    epsilon = 10 ** -12
    X = np.linspace(0 + epsilon, 1, samples)
    X_eval = [pysgpp.DataVector([i]) for i in X]

    return X, X_eval


def calculate_grid_y_values(X, operation, level):
    X_eval = [pysgpp.DataVector([i]) for i in X]
    Y = [operation.evaluate(level, X_eval[i]) for i in range(len(X_eval))]

    z = operation.getLevelManager().getGridPointMatrix()

    gridpoints = [z.get(0, x) for x in range(z.getNcols())]

    return Y, gridpoints


def calc_error_ilevels(operation, targetfunc, dim, maxlevel):
    nr_gridpoints = []
    errors = []
    for l in range(1, maxlevel):
        operation_wrap = operationwrapper(operation, l)
        gridpterror = estimatel2Error(10000, dim, operation_wrap, targetfunc)
        errors.append(gridpterror)
        nr_gridpoints.append(operation.getLevelManager().numGridPoints())
    return nr_gridpoints, errors


def calc_error_ilevels_grad(operation, func_container, dim, maxlevel, grad_index_list):
    nr_gridpoints = []
    errors = []

    for l in range(1, maxlevel):
        operation_wrap = operationwrapper(operation, l)
        for k in grad_index_list:
            operation_wrap = fctClass.getgradkfunc(operation_wrap, k)
        targetfunc = func_container.getGradient(grad_index_list)

        gridpterror = estimatel2Error(10000, dim, operation_wrap, targetfunc)
        errors.append(gridpterror)
        nr_gridpoints.append(operation.getLevelManager().numGridPoints())
    return nr_gridpoints, errors


def plot_log_error(errors, gridpoints, text, x_axis_type="nr_gridpoints", title="", filename=None,
                   show=True):
    fig = plt.figure()
    plt.title(title)
    ax = fig.add_subplot(111)
    for i in range(errors.shape[1]):
        ax.plot(gridpoints[:, i], errors[:, i], marker='o', label=text[i])

    plt.legend()
    if (x_axis_type == "nr_gridpoints"):
        ax.set_xlabel("Anzahl Gitterpunkte")
    else:
        ax.set_xlabel("Level")
    ax.set_ylabel("$L^2-Fehler$")
    if (x_axis_type == "nr_gridpoints"):
        ax.set_xscale("log")
    ax.set_yscale("log")

    plt.grid(b=True, which='major', color='k', linestyle='--')

    if (show == True):
        plt.show()

    if filename != None:
        plt.savefig(filename + '.pdf', format='pdf', dpi=900)
        plt.close("all")


def example_1D_realfunction():
    function = f1D
    X, X_eval = generate1DGrid(500)
    Y_compare = [function(X_eval[i]) for i in range(len(X_eval))]
    p.plot_1d_grid(X, Y_compare, "function", )


def example_1D_psi(l):
    func = pysgpp.multiFunc(f1D)
    d = 1
    level = l
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformnakBsplineInterpolation(
        d, func)
    X, _, = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    p.plot_1d_grid(X, Y, "$\psi$", gridpoints)


def example_1D_linear(l):
    func = pysgpp.multiFunc(f1D)
    d = 1
    level = 1
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        d, func)
    X, _, = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    p.plot_1d_grid(X, Y, "$\linear$", gridpoints)


def example_1D_zeta(l):
    func = pysgpp.multiFunc(fctClass.getgradkfunc(f1D, 0))
    d = 1
    level = l
    n_samples = 500
    operation = pysgpp.CombigridOperation.createExpUniformBoundaryZetaHermiteInterpolation(
        d, func)
    X, _, = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    p.plot_1d_grid(X, Y, "$\psi$", gridpoints)


def example_combicombigrid_1D(l):
    x = pysgpp.DataVector([0.01])
    x1 = pysgpp.DataVector([0.5])
    x2 = pysgpp.DataVector([0.99])

    d = 1
    level = l

    f1D_grad = fctClass.getgradkfunc(f1D, 0)
    f1D_grad_wrap = pysgpp.multiFunc(f1D_grad)
    f1D_wrap = pysgpp.multiFunc(f1D)
    operation = gC.CombiCombigrid1d_hermite(f1D, f1D_grad)
    op_wrap = operationwrapper(operation, level)

    grad_operation = fctClass.getgradkfunc(op_wrap, 0)

    print(str(x), grad_operation(x), "ableitung")
    print(str(x1), grad_operation(x1), "ableitung")
    print(str(x1), f1D_grad(x1), "ableitung real")
    print(str(x2), grad_operation(x2), "ableitung")

    # plot the grid
    X, X_eval = generate1DGrid(500)
    Y, gridpoints = calculate_grid_y_values(X, operation, level)
    print(
        str(operation.operation_zeta.numGridPoints()) + " grid points in zeta/psi")
    p.plot_1d_grid(X, Y, "combi_combigrid_1D", gridpoints)

    # error calculation

    error = estimatel2Error(
        7000, 1, operationwrapper(operation, level), f1D)
    print("the error is " + str(error))


def example_2D_linear(l, func_standard):
    d = 2
    level = l
    n_samples = 100

    operation = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        d, func_standard)

    operation_wrap = operationwrapper(operation, level)
    p.plot2DGrid_with_tangents(n_samples, level, operation, "linear")

    print(operation.evaluate(level, pysgpp.DataVector([0, 1])))
    print(operation.evaluate(level, pysgpp.DataVector([0.000000001, 1])))

    error = estimatel2Error(20000, 2, operation_wrap, func_standard)
    # calculate error
    print("the estimtaded l2 error for lineargrid: " + str(error))


def example_2D_psi():
    d = 2
    level = 2
    n_samples = 100
    func = pysgpp.OptRosenbrockObjective(d)
    func_standard = getfuncwrapper(func)
    func_combiwrap = pysgpp.multiFunc(func_standard)

    operation = pysgpp.CombigridOperation.createExpUniformBoundaryPsiHermiteInterpolation(
        d, func_combiwrap)

    operation_wrap = operationwrapper(operation, level)
    p.plot2DGrid_operation(n_samples, level, operation, func_combiwrap)

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)
    # calculate error
    print("the estimtaded l2 error for psigrid: " + str(error))


def example_2D_comparison_function(func_standard, title, show=True):
    d = 2

    n_samples = 100

    X = np.linspace(0, 1, n_samples)
    Y = np.linspace(0, 1, n_samples)

    p.plot2D_comparison_function(n_samples, func_standard, title, show=show)


def example_combicombigrid_2D_linear(l, func_standard):
    """

    Args:
        l:
        func_standard: function that is already in multifunc standard format
    """
    d = 2
    level = 0
    n_samples = 50

    operation = pysgpp.CombigridOperation.createExpUniformnakBsplineInterpolation(
      d, func_standard)
    operation_wrap = operationwrapper(operation, level)

    p.plot2DGrid_operation(n_samples, level, operation, "de Baar & Harding")
    p.plot2DContour(n_samples, level, operation)

    # derivatives


    grad_0 = fctClass.getgradkfunc(operation_wrap, 0)
    grad_1 = fctClass.getgradkfunc(operation_wrap, 1)

    x1 = pysgpp.DataVector([0.5, 0.5])
    # print(grad_0(x1))
    # print(grad_1(x1))


    error = estimatel2Error(20000, 2, operation_wrap, func_standard)

    # calculate error
    print("the estimtaded l2 error for combicombilinear: " + str(error))

    gridpterror = calc_error_on_gridpts(operation.getLevelManager().getAllGridPoints(),
                                        operation_wrap, func_standard)
    grad_error = calc_error_gradient(operation.getLevelManager().getAllGridPoints(), operation_wrap,
                                     func_standard, d)
    mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
                                              operation_wrap,
                                              func_standard, d)

    print("estimtaded l2 error for combicombihermite: " + str(error))
    print("1st gradient errors:" + str(grad_error))
    print("mixed gradient error:" + str(mixgrad_error))


def example_combicombigrid_2D(l, func_collection):
    d = 2
    level =2
    n_samples = 50
    func_standard = func_collection.getFunction()
    #operation = gC.CombiCombigriddeBaarHarding(func_collection, 2)
    operation=pysgpp.CombigridOperation.createExpUniformBoundaryBsplinePsiInterpolation(d,0,\
              pysgpp.multiFunc(func_standard),3)
    # operation = HierachGridBSpline(2, 3, func_standard)
    operation_wrap = operationwrapper(operation, level)

    # operation.operation_zeta for mixed
    p.plot2DGrid_operation(n_samples, level, operation, "Psi-BSPline",show=False)
    #p.plot2DContour(n_samples, level, operation)
    plt.show()

    error = estimatel2Error(10000, 2, operation_wrap, func_standard)

    # grad_error = calc_error_gradient(operation.getLevelManager().getAllGridPoints(),
    #                                 operation_wrap,
    #                                 func_standard, d)
    # mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
    #                                          calc_error_mixed_gradient,
    #                                          func_standard, d)
    # calculate error
    print("estimtaded l2 error for combicombihermite: " + str(error))
    # print("1st gradient errors:" + str(grad_error))
    # print("mixed gradient error:" + str(mixgrad_error))


def geterrorfunc(func1, func2):
    def errorfunc(x):
        return abs(func1(x) - func2(x))

    return errorfunc


def example_error_picewise(l, func_container, show=True):
    d = 2
    level = l
    n_samples = 50

    operation_list = []

    func_standard = func_container.getFunction()
    func_standard = pysgpp.multiFunc(func_standard)

    operation_list.append(gC.CombiCombigriddeBaarHarding(func_container, 2))
    operation_list.append(gC.CombiCombigridHermite(func_container, 2))
    operation_list.append(pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        2, func_standard))
    operation_list.append(gC.CombiFullGridHermite(func_container, 2))
    operation_list.append(gC.CombiFullGridLinear(func_container, 2))
    # operation_list.append(BaseLinearFullgrid(2, func_container, True))
    operation_list.append(gC.HierachGridBSpline(2, 3, func_container))

    text = ["Combilinear", "combiHermite", "linear", "CombiHermiteFull", "CombiLinearFull",
            "BSpline(3)"
            ]

    i = 0
    for operation in operation_list:
        operation_wrap = operationwrapper(operation, level)
        operation_container = fctClass.funcGradientCollection(operation_wrap, dim)

        # operation.operation_zeta for mixed
        # plot2D_comparison_function_random(n_samples, geterrorfunc(operation_wrap, func_standard),
        #                                 text[i])


        p.plot2DContour_func(n_samples, geterrorfunc(operation_wrap, func_standard),
                             operation, title=text[i])

        error = estimatel2Error(10000, 2, operation_wrap, func_standard)
        errorl2_gradients = estimatel2ErrorGradients(10000, 2, operation_container, func_container)

        erroron_gripoint = calc_error_on_gridpts(operation.getLevelManager().getAllGridPoints(),
                                                 operation_wrap, func_standard)
        mixgrad_error = calc_error_mixed_gradient(operation.getLevelManager().getAllGridPoints(),
                                                  operation_container,
                                                  func_container, d)

        print("estimtaded l2 error for " + text[i] + ": " + str(error))
        print("error on gridpts: " + str(erroron_gripoint))
        print("mixed gradient error:" + str(mixgrad_error))
        print("l2 errors from gradients: " + str(errorl2_gradients))

        print("-----------------------------------------------------------------")
        i = i + 1
    if (show == True):
        plt.show()


def example_plot_error(func_collection, dim, maxlevel=5, x_axis_type="nr_gridpoints", title="",
                       show=True,
                       filename=None):

    """

      Args:
          func_standard:
          dim:
          maxlevel:
          x_axis_type: either "nr_gridpoints" for the number of gridpoints or "level" to compare
          each operation by level
      """

    # operation4= HierachGridBSpline(dim,3,func_standard)


    func_standard = pysgpp.multiFunc(func_collection.getFunction())

    operations = []

    operations.append(gC.CombiCombigriddeBaarHarding(func_collection, dim))
    operations.append(gC.CombiCombigridHermite(func_collection, dim))
    operations.append(pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        dim, func_standard))
    operations.append(gC.CombiCombigrid2dHermite_without_mixed(func_collection, dim))
    # operations.append((BaseLinearFullgrid(dim, func_collection, True)))

    # operations.append(gC.LinearFullgrid(func_collection, dim))
    # operations.append(CombiFullGridHermite(func_collection, dim))
    # operations.append(CombiFullGridLinear(func_collection, dim))
    # operations.append(HierachGridBSpline(dim, 3, func_collection, True))
    operations.append(gC.HierachGridBSpline(dim, 3, func_collection, False))
    operations.append(pysgpp.CombigridOperation.createExpUniformBoundaryBsplineInterpolation(
        dim, func_standard, 3))

    func_standard = func_collection.getFunction()
    labels = []
    labels.append("Gitter (de Baar und Harding)")
    labels.append("Gitter (Hermite)")
    labels.append("lineares Gitter")
    labels.append("Gitter (Hermite) ohne gemischte Abl.")

    # labels.append("BaseLinearFullgrid")
    # labels.append("linearFullgrid")
    # labels.append("HermiteFullgrid")
    # labels.append("CombiLinearFullgrid")
    # labels.append("BSpline-Grid-Full")
    labels.append("Base_BSpline-Grid(Grad 3)")
    labels.append("BSpline-Grid(Grad 3)")

    X = np.zeros((maxlevel - 1, len(operations)))
    Y = np.zeros((maxlevel - 1, len(operations)), dtype=np.float64)

    for i in range(len(operations)):
        nr_gridpoints, errors = calc_error_ilevels(operations[i], func_standard, dim, maxlevel)

        if (x_axis_type == "level"):
            x_values = [j for j in range(1, maxlevel)]

        else:
            x_values = nr_gridpoints
        print(i, nr_gridpoints)

        X[:, i] = x_values
        Y[:, i] = errors

        if (filename is not ""):

            sl.save_namelist(labels, "labels.txt", filename)
            sl.save_array(X, "X_values.data", filename)
            sl.save_array(Y, "Y_values.data", filename)

        else:
            print("no filename")


def loadandplot(X_name, Y_name, dirname, labelname, title, x_axis_type="nr_gridpoints", show=True,
                filename=None):
    X = sl.load_array(X_name, dirname)
    Y = sl.load_array(Y_name, dirname)
    labels = sl.load_list(labelname, dirname)

    plot_log_error(Y, X, labels, x_axis_type, title=title, filename=filename, show=show)


def example_plot_error_gradients(func_collection, dim, grad_index_list, maxlevel=5,
                                 x_axis_type="nr_gridpoints",
                                 title="", filename=None, show=True):


    """

      Args:
          func_standard:
          dim:
          maxlevel:
          x_axis_type: either "nr_gridpoints" for the number of gridpoints or "level" to compare
          each operation by level
      """




    func_standard = pysgpp.multiFunc(func_collection.getFunction())

    operations = []

    operations.append(gC.CombiCombigriddeBaarHarding(func_collection, dim))
    operations.append(gC.CombiCombigridHermite(func_collection, dim))
    operations.append(pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(
        dim, func_standard))
    operations.append(gC.CombiCombigrid2dHermite_without_mixed(func_collection, dim))
    # operations.append((BaseLinearFullgrid(dim, func_collection, True)))

    # operations.append(gC.LinearFullgrid(func_collection, dim))
    # operations.append(gC.CombiFullGridHermite(func_collection, dim))
    # operations.append(CombiFullGridLinear(func_collection, dim))
    # operations.append(gC.CombiFullGridHermite_withoutmixed(func_collection, dim))
    # operations.append(HierachGridBSpline(dim, 3, func_collection, True))
    operations.append(gC.HierachGridBSpline(dim, 3, func_collection, False))
    operations.append(pysgpp.CombigridOperation.createExpUniformBoundaryBsplineInterpolation(
          dim,func_standard,3))


    labels = []
    labels.append("Gitter (de Baar und Harding)")
    labels.append("Gitter (Hermite)")
    labels.append("lineares Gitter")
    labels.append("Gitter (Hermite) ohne gemischte Abl.")

    # labels.append("BaseLinearFullgrid")
    # labels.append("linearFullgrid")
    # labels.append("HermiteFullgrid")
    # labels.append("CombiLinearFullgrid")
    # labels.append("BSpline-Grid-Full")
    labels.append("Base_BSpline-Grid(Grad 3)")
    labels.append("BSpline-Grid(Grad 3)")

    X = np.zeros((maxlevel - 1, len(operations)))
    Y = np.zeros((maxlevel - 1, len(operations)), dtype=np.float64)

    for i in range(len(operations)):
        nr_gridpoints, errors = calc_error_ilevels_grad(operations[i], func_collection, dim,
                                                        maxlevel, grad_index_list)

        if (x_axis_type == "level"):
            x_values = [j for j in range(1, maxlevel)]

        else:
            x_values = nr_gridpoints
        print(i, nr_gridpoints)

        X[:, i] = x_values
        Y[:, i] = errors

        if (filename is not ""):

            sl.save_namelist(labels, "labels.txt", filename)
            sl.save_array(X, "X_values.data", filename)
            sl.save_array(Y, "Y_values.data", filename)

        else:
            print("no filename")
            # plot_log_error(Y, X, labels, x_axis_type, title=title, filename=filename, show=show)


def example_calcl2error(func_container, name, level, dim, fct=True, grad_index=[]):
    filename = name
    if (fct):
        example_plot_error(func_container, dim, maxlevel=level, x_axis_type="nr_gridpoints",
                           title=filename, filename=filename,
                           show=False)

    for i in grad_index:
        filename = name

        for j in i:
            filename += "x" + str(j)
        example_plot_error_gradients(func_container, dim, i, maxlevel=level,
                                     x_axis_type="nr_gridpoints",
                                     title=filename, filename=filename,
                                     show=False)


def plot_l2error(name, fct=True, grad_index=[],x_axis_type="nr_gridpoints"):
    filename = name
    if (fct):
        loadandplot("X_values.data", "Y_values.data", filename, "labels.txt", title=filename,
                    x_axis_type=x_axis_type
                    , show=True)

    for i in grad_index:
        filename = name

        for j in i:
            filename += "x" + str(j)

        loadandplot("X_values.data", "Y_values.data", filename, "labels.txt", title=filename,
                    x_axis_type=x_axis_type
                    , show=True)


def testf(x):
    return x[0] * x[1]


def examples_1d():
    # example_1D_realfunction()
    # example_1D_psi(1)

    # example_1D_linear(1)

    # example_1D_zeta(2)
    # example_combicombigrid_1D(2)
    plt.show()

def examples_2d():

    dim = 2

    func = pysgpp.OptBraninObjective()
    func_wrap = getfuncwrapper(func)

    func_standard = pysgpp.multiFunc(func_wrap)

    # example_2D_psi()
    # example_2D_linear(2,func_standard)
    #example_combicombigrid_2D_linear(2, func_standard)  # with "contourplot"
    # print()

    testclass = fctClass.funcGradientCollection(func_standard, 2)

    func_container = fctClass.funcGradientCollectionSymbolic(fctClass.BraninSymbolic(), 2)

    #example_2D_comparison_function(func_container.getFunction(), "", show=False)

    #example_combicombigrid_2D(3, func_container)

    #example_calcl2error(func_container, "Values/Test/Branin/BSplines", 9, 2, grad_index=[[0], [1],
    #                                                                                   [0, 1]],
    #                  fct=True)
    plot_l2error("Values/Test/Branin/BSplines",grad_index=[[0],[1],[0,1]],fct=True)




    #example_error_picewise(2, func_container)


def examples_multid():
    dim = 4

    func_container = fctClass.funcGradientCollectionSymbolic(fctClass.testfSymbolic2_4d(), dim)


    # plot_l2error("Values/test4d",grad_index=[[1],[1],[0,1]],fct=True)
    # example_calcl2error(func_container, "test4d", level=2,grad_index=[[0, 1], [1]],fct=True)

    # func = pysgpp.OptRosenbrockObjective(dim)
    # func_wrap = getfuncwrapper(func)

    # func_standard = pysgpp.multiFunc(func_wrap)


#examples_1d()

examples_2d()

#examples_multid()




// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_constrainedOptimization_cpp constrainedOptimization.cpp
 * This example demonstrates the optimization of an objective function \f$ f\f$ with additional
 * constraints.
 * The inequality constraints are specified via a function \f$ g\f$, the equality constraints are
 * specified via a function \f$ h\f$.
 */
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <sgpp/optimization/test_problems/constrained/Floudas.hpp>
#include <sgpp/optimization/test_problems/constrained/G03.hpp>
#include <sgpp/optimization/test_problems/constrained/G04.hpp>
#include <sgpp/optimization/test_problems/constrained/G05.hpp>
#include <sgpp/optimization/test_problems/constrained/G06.hpp>
#include <sgpp/optimization/test_problems/constrained/G08.hpp>
#include <sgpp/optimization/test_problems/constrained/G09.hpp>
#include <sgpp/optimization/test_problems/constrained/G10.hpp>
#include <sgpp/optimization/test_problems/constrained/G11.hpp>
#include <sgpp/optimization/test_problems/constrained/G12.hpp>
#include <sgpp/optimization/test_problems/constrained/G13.hpp>
#include <sgpp/optimization/test_problems/constrained/Simionescu.hpp>
#include <sgpp/optimization/test_problems/constrained/Soland.hpp>

#include <iostream>

/**
 * Set the verbosity of the printed output. Increase the argument for more information about the
 * process of solving linear system.
 */
int main() {
  sgpp::base::Printer::getInstance().setVerbosity(0);

  /**
   * Load the predefined test problem G04. G04 consists of the following:
   * - Number of parameters: 5\n
   * - Number of inequality constraints: 6\n
   * - Number of equality constraints: 0\n
   * - Domain: \f$\bar{\vec{x}} \in  [78, 102] \times [33, 45] \times [27, 45]^3\f$\n
   * - Optimal point: \f$\bar{\vec{x}}_{\text{opt}} = (78, 33, 29.99526, 45, 36.77581)\f$\n
   * - Optimal function value: \f$\bar{f}(\bar{\vec{x}}_{\text{opt}}) =  -30665.54\f$\n
   */
  sgpp::optimization::test_problems::G04 problem;

  problem.generateDisplacement();

  sgpp::base::ScalarFunction& f = problem.getObjectiveFunction();
  sgpp::base::VectorFunction& g = problem.getInequalityConstraintFunction();
  sgpp::base::VectorFunction& h = problem.getEqualityConstraintFunction();

  /**
   * Generate a regular sparse grid of level five with modified B-Splines of degree five as basis
   * functions.
   */
  const size_t p = 5;
  const size_t d = problem.getObjectiveFunction().getNumberOfParameters();
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createModBsplineGrid(d, p));
  const size_t l = 5;
  grid->getGenerator().regular(l);

  /**
   * Initialize vector \f$ x \f$ that will contain the grid points and auxiliary vectors
   * that will contain \f$ g(x) \f$ and \f$ h(x) \f$ for the hierarchisation.  Initialize further a
   * vector and two matrices that will contain the hierarchisation coefficients \f$ f_{\alpha} ,
   * g_{\alpha} \f$ and \f$ h_{\alpha} \f$.
   */
  const size_t N = grid->getSize();
  const size_t mG = g.getNumberOfComponents();
  const size_t mH = h.getNumberOfComponents();
  sgpp::base::DataVector x(d), gx(mG), hx(mH);
  sgpp::base::DataVector fAlpha(N);
  sgpp::base::DataMatrix gAlpha(N, mG);
  sgpp::base::DataMatrix hAlpha(N, mH);

  /**
   * Prepare the hierarchisation by filling  \f$ f_{\alpha}, g_{\alpha}\f$ and \f$ h_{\alpha} \f$
   * with the function values \f$ f(x), g(x)\f$ and \f$ h(x)\f$ for every grid point x.
   */
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  for (size_t i = 0; i < N; i++) {
    x = gridStorage.getCoordinates(gridStorage[i]);
    fAlpha[i] = f.eval(x);
    g.eval(x, gx);
    gAlpha.setRow(i, gx);
    h.eval(x, hx);
    hAlpha.setRow(i, hx);
  }

  /**
   * Perform the hierarchisation.
   */
  std::unique_ptr<sgpp::optimization::OperationMultipleHierarchisation> hierOp(
      sgpp::op_factory::createOperationMultipleHierarchisation(*grid));
  hierOp->doHierarchisation(fAlpha);
  hierOp->doHierarchisation(gAlpha);
  hierOp->doHierarchisation(hAlpha);

  /**
   * Using the coefficients \f$ f_{\alpha}, g_{\alpha}\f$ and \f$ h_{\alpha} \f$ define the
   * interpolant functions \f$\tilde{f}, \tilde{g}\f$ and \f$\tilde{h}\f$ as well as their
   * gradients.
   */
  sgpp::base::InterpolantScalarFunction ft(*grid, fAlpha);
  sgpp::base::InterpolantScalarFunctionGradient ftGradient(*grid, fAlpha);

  sgpp::base::InterpolantVectorFunction gt(*grid, gAlpha);
  sgpp::base::InterpolantVectorFunctionGradient gtGradient(*grid, gAlpha);

  sgpp::base::InterpolantVectorFunction ht(*grid, hAlpha);
  sgpp::base::InterpolantVectorFunctionGradient htGradient(*grid, hAlpha);

  /**
   * Define an augmented lagrangian method to solve the constrained optimization problem using the
   * interpolants \f$\tilde{f}, \tilde{g}\f$ and \f$\tilde{h}\f$ as well as their gradients as
   * arguments.
   */
  sgpp::optimization::optimizer::AugmentedLagrangian optimizer(ft, ftGradient, gt, gtGradient, ht,
                                                               htGradient);

  optimizer.setN(10000);

  /**
   * Create a feasible starting point \f$x_0 \f$ and apply the augmented lagrangian method defined
   * above.
   */
  const sgpp::base::DataVector x0(optimizer.findFeasiblePoint());
  optimizer.setStartingPoint(x0);
  optimizer.optimize();

  /**
   * Get the optimum point \f$ x_{Opt}\f$ and the calculated approximation to the optimum point \f$
   * x_{Opt}^*\f$.
   */
  sgpp::base::DataVector xOpt(d);
  problem.getOptimalPoint(xOpt);
  const sgpp::base::DataVector xOptStar(optimizer.getOptimalPoint());

  /**
   * Evaluate \f$ f,g\f$ and \f$h\f$ as well as their gradients in \f$ \lbrace x, x_{Opt},
   * x_{Opt}^{*} \rbrace \f$
   */
  for (const sgpp::base::DataVector& x : {xOpt, x0, xOptStar}) {
    const double fx = f.eval(x);
    sgpp::base::DataVector gx(mG);
    g.eval(x, gx);
    sgpp::base::DataVector hx(mH);
    h.eval(x, hx);

    const double ftx = ft.eval(x);
    sgpp::base::DataVector gtx(mG);
    gt.eval(x, gtx);
    sgpp::base::DataVector htx(mH);
    ht.eval(x, htx);

    /**
     * Print the evaluation point, and the values of \f$ f, g, h\f$ as well as the values
     * \f$\tilde{f}, \tilde{g}\f$ and \f$\tilde{h}\f$ evaluated in the evaluation point.
     */
    std::cout << "x = " << x.toString() << "\n";
    std::cout << "f = " << fx << ", g = " << gx.toString() << ", h = " << hx.toString() << "\n";
    std::cout << "ft = " << ftx << ", gt = " << gtx.toString() << ", ht = " << htx.toString()
              << "\n";
    std::cout << "\n";
  }

  return 0;
}
/**
 * When executed this example produces the following output\n
 *  (remark that the used example problem G04 has no equality constraints, so h is empty) \n
 *  \n
 *  \n
 * Solving linear system (automatic method)...\n
Done in 11124ms.\n
Solving linear system (automatic method)...\n
Done in 4013ms.\n
Solving linear system (automatic method)...\n
Done in 3746ms.\n
Optimizing (Augmented Lagrangian)...\n
Done in 157ms.\n
Optimizing (Augmented Lagrangian)...\n
Done in 7290ms.\n
x = [8.07474373435067976912e-03, 6.41782974255440154254e-03, 1.61029433421242457181e-01,
9.96927758912899086852e-01, 5.50146141071411132195e-01]\n
f = -30665.5, g = [-4.26325641456060111523e-14, -9.19999999999999573674e+01,
-1.11594996910731083517e+01, -8.84050030892689164830e+00, -4.99999999999988631316e+00,
-1.13686837721616029739e-13], h = []\n
ft = -30665.5, gt = [2.97498166580675220360e-14, -9.20000000000000710543e+01,
-1.11594988391265399486e+01, -8.84050116087346715688e+00, -4.99999999999987476684e+00,
-1.14863590187006103699e-13], ht = []\n
\n
x = [4.19202291911377755707e-01, 1.35037410213593056518e-01, 7.12332125450077957574e-01,
3.49062366127655510084e-01, 1.19669759364614158859e-01]\n
f = -26846.1, g = [-1.68542663525514058165e+00, -9.03145733647448594184e+01,
-9.76861095224427344874e+00, -1.02313890477557265513e+01, -3.30920437708866188586e+00,
-1.69079562291133811414e+00], h = []\n
ft = -26846.1, gt = [-1.68542663525515057366e+00, -9.03145733647449446835e+01,
-9.76861091582623330964e+00, -1.02313890841737737958e+01, -3.30920437708865167181e+00,
-1.69079562291134677388e+00], ht = []\n
\n
x = [0.00000000000000000000e+00, 1.33215347999698752179e-01, 6.77166052213670210946e-01,
9.99970269638240760735e-01, 9.02152104816860589409e-01]\n
f = -26818.1, g = [2.56624090454906195191e-01, -9.22566240904549061952e+01,
-7.46102616932890327917e+00, -1.25389738306710967208e+01, -5.21058637424896886614e-01,
-4.47894136257510311339e+00], h = []\n
ft = -26818.1, gt = [2.56624090454927789029e-01, -9.22566240904548919843e+01,
-7.46102618953152685322e+00, -1.25389738104684500541e+01, -5.21058637424886450518e-01,
-4.47894136257504182907e+00], ht = []

 *
 *
 */

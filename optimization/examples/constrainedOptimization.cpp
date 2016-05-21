// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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

int main() {
  std::cout << "Hello World!\n";

  sgpp::optimization::Printer::getInstance().setVerbosity(5);
  sgpp::optimization::RandomNumberGenerator::getInstance().setSeed();

  sgpp::optimization::test_problems::G04 problem;
  const size_t d = problem.getObjectiveFunction().getNumberOfParameters();

  /*const size_t d = 6;
  sgpp::optimization::test_problems::ProductOnSphere problem(d);*/

  problem.generateDisplacement();

  sgpp::optimization::ScalarFunction& f =
    problem.getObjectiveFunction();
  sgpp::optimization::VectorFunction& g =
    problem.getInequalityConstraintFunction();
  sgpp::optimization::VectorFunction& h =
    problem.getEqualityConstraintFunction();

  /*sgpp::optimization::WrapperScalarFunctionGradient fGradient(
    d, [d](const sgpp::base::DataVector & x,
  sgpp::base::DataVector & gradient) {
    double fx = -1.0;
    gradient.setAll(-1.0);

    for (size_t t = 0; t < d; t++) {
      for (size_t t2 = 0; t2 < d; t2++) {
        if (t2 != t) {
          gradient[t2] *= x[t];
        } else {
          fx *= x[t];
        }
      }
    }

    return fx;
  });

  sgpp::optimization::WrapperVectorFunctionGradient gGradient(
    d, 1, [d](const sgpp::base::DataVector & x,
              sgpp::base::DataVector & value,
  sgpp::base::DataMatrix & gradient) {
    double gx = -1.0;

    for (size_t t = 0; t < d; t++) {
      gx += x[t] * x[t];
      gradient(0, t) = 2.0 * x[t];
    }

    value[0] = gx;
  });

  sgpp::optimization::VectorFunctionGradient& hGradient =
    sgpp::optimization::EmptyVectorFunctionGradient::getInstance();*/

  const size_t mG = g.getNumberOfComponents();
  const size_t mH = h.getNumberOfComponents();

  const size_t p = 5;
  std::unique_ptr<sgpp::base::Grid> grid(
    sgpp::base::Grid::createModBsplineGrid(d, p));

  const size_t l = 5;
  grid->getGenerator().regular(l);

  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  const size_t N = grid->getSize();

  std::unique_ptr<sgpp::optimization::OperationMultipleHierarchisation> hierOp(
    sgpp::op_factory::createOperationMultipleHierarchisation(*grid));

  sgpp::base::DataVector x(d), gx(mG), hx(mH);
  sgpp::base::DataVector fAlpha(N);
  sgpp::base::DataMatrix gAlpha(N, mG);
  sgpp::base::DataMatrix hAlpha(N, mH);

  for (size_t i = 0; i < N; i++) {
    const sgpp::base::GridIndex& gp = gridStorage[i];

    for (size_t t = 0; t < d; t++) {
      x[t] = gp.getCoord(t);
    }

    fAlpha[i] = f.eval(x);

    g.eval(x, gx);
    gAlpha.setRow(i, gx);

    h.eval(x, hx);
    hAlpha.setRow(i, hx);
  }

  hierOp->doHierarchisation(fAlpha);
  hierOp->doHierarchisation(gAlpha);
  hierOp->doHierarchisation(hAlpha);

  sgpp::optimization::InterpolantScalarFunction ft(*grid, fAlpha);
  sgpp::optimization::InterpolantScalarFunctionGradient ftGradient(
    *grid, fAlpha);

  sgpp::optimization::InterpolantVectorFunction gt(*grid, gAlpha);
  sgpp::optimization::InterpolantVectorFunctionGradient gtGradient(
    *grid, gAlpha);

  sgpp::optimization::InterpolantVectorFunction ht(*grid, hAlpha);
  sgpp::optimization::InterpolantVectorFunctionGradient htGradient(
    *grid, hAlpha);

  sgpp::optimization::optimizer::AugmentedLagrangian optimizer(
    ft, ftGradient, gt, gtGradient, ht, htGradient);
  /*sgpp::optimization::optimizer::AugmentedLagrangian optimizer(
    f, fGradient, g, gGradient, h, hGradient);*/

  optimizer.setN(10000);

  const sgpp::base::DataVector x0(optimizer.findFeasiblePoint());
  // const sgpp::base::DataVector x0(d, 0.2);

  optimizer.setStartingPoint(x0);
  optimizer.optimize();

  sgpp::base::DataVector xOpt(d);
  problem.getOptimalPoint(xOpt);
  const sgpp::base::DataVector xOptStar(optimizer.getOptimalPoint());

  for (const sgpp::base::DataVector& x : {
  xOpt, x0, xOptStar
}) {
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

    std::cout << "x = " << x.toString() << "\n";
    std::cout << "f = " << fx << ", g = " << gx.toString() << ", h = " <<
              hx.toString() << "\n";
    std::cout << "ft = " << ftx << ", gt = " << gtx.toString() << ", ht = " <<
              htx.toString() << "\n";
    std::cout << "\n";
  }

  return 0;
}

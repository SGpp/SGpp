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

  SGPP::optimization::Printer::getInstance().setVerbosity(5);
  SGPP::optimization::RandomNumberGenerator::getInstance().setSeed();

  SGPP::optimization::test_problems::G04 problem;
  const size_t d = problem.getObjectiveFunction().getNumberOfParameters();

  /*const size_t d = 6;
  SGPP::optimization::test_problems::ProductOnSphere problem(d);*/

  problem.generateDisplacement();

  SGPP::optimization::ScalarFunction& f =
    problem.getObjectiveFunction();
  SGPP::optimization::VectorFunction& g =
    problem.getInequalityConstraintFunction();
  SGPP::optimization::VectorFunction& h =
    problem.getEqualityConstraintFunction();

  /*SGPP::optimization::WrapperScalarFunctionGradient fGradient(
    d, [d](const SGPP::base::DataVector & x,
  SGPP::base::DataVector & gradient) {
    SGPP::float_t fx = -1.0;
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

  SGPP::optimization::WrapperVectorFunctionGradient gGradient(
    d, 1, [d](const SGPP::base::DataVector & x,
              SGPP::base::DataVector & value,
  SGPP::base::DataMatrix & gradient) {
    SGPP::float_t gx = -1.0;

    for (size_t t = 0; t < d; t++) {
      gx += x[t] * x[t];
      gradient(0, t) = 2.0 * x[t];
    }

    value[0] = gx;
  });

  SGPP::optimization::VectorFunctionGradient& hGradient =
    SGPP::optimization::EmptyVectorFunctionGradient::getInstance();*/

  const size_t mG = g.getNumberOfComponents();
  const size_t mH = h.getNumberOfComponents();

  const size_t p = 5;
  std::unique_ptr<SGPP::base::Grid> grid(
    SGPP::base::Grid::createModBsplineGrid(d, p));

  const size_t l = 5;
  std::unique_ptr<SGPP::base::GridGenerator>(
    grid->createGridGenerator())->regular(l);

  SGPP::base::GridStorage& gridStorage = grid->getStorage();

  const size_t N = grid->getSize();

  std::unique_ptr<SGPP::optimization::OperationMultipleHierarchisation> hierOp(
    SGPP::op_factory::createOperationMultipleHierarchisation(*grid));

  SGPP::base::DataVector x(d), gx(mG), hx(mH);
  SGPP::base::DataVector fAlpha(N);
  SGPP::base::DataMatrix gAlpha(N, mG);
  SGPP::base::DataMatrix hAlpha(N, mH);

  for (size_t i = 0; i < N; i++) {
    const SGPP::base::GridIndex& gp = *gridStorage[i];

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

  SGPP::optimization::InterpolantScalarFunction ft(*grid, fAlpha);
  SGPP::optimization::InterpolantScalarFunctionGradient ftGradient(
    *grid, fAlpha);

  SGPP::optimization::InterpolantVectorFunction gt(*grid, gAlpha);
  SGPP::optimization::InterpolantVectorFunctionGradient gtGradient(
    *grid, gAlpha);

  SGPP::optimization::InterpolantVectorFunction ht(*grid, hAlpha);
  SGPP::optimization::InterpolantVectorFunctionGradient htGradient(
    *grid, hAlpha);

  SGPP::optimization::optimizer::AugmentedLagrangian optimizer(
    ft, ftGradient, gt, gtGradient, ht, htGradient);
  /*SGPP::optimization::optimizer::AugmentedLagrangian optimizer(
    f, fGradient, g, gGradient, h, hGradient);*/

  optimizer.setN(10000);

  const SGPP::base::DataVector x0(optimizer.findFeasiblePoint());
  // const SGPP::base::DataVector x0(d, 0.2);

  optimizer.setStartingPoint(x0);
  optimizer.optimize();

  SGPP::base::DataVector xOpt(d);
  problem.getOptimalPoint(xOpt);
  const SGPP::base::DataVector xOptStar(optimizer.getOptimalPoint());

  for (const SGPP::base::DataVector& x : {
  xOpt, x0, xOptStar
}) {
    const SGPP::float_t fx = f.eval(x);
    SGPP::base::DataVector gx(mG);
    g.eval(x, gx);
    SGPP::base::DataVector hx(mH);
    h.eval(x, hx);

    const SGPP::float_t ftx = ft.eval(x);
    SGPP::base::DataVector gtx(mG);
    gt.eval(x, gtx);
    SGPP::base::DataVector htx(mH);
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

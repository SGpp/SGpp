// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/optimization/test_problems/unconstrained/AbsoluteValue.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Ackley.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Alpine02.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Beale.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Branin01.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Branin02.hpp>
#include <sgpp/optimization/test_problems/unconstrained/BubbleWrap.hpp>
#include <sgpp/optimization/test_problems/unconstrained/EasomYang.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Eggholder.hpp>
#include <sgpp/optimization/test_problems/unconstrained/GoldsteinPrice.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Griewank.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Hartman3.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Hartman6.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Himmelblau.hpp>
#include <sgpp/optimization/test_problems/unconstrained/HoelderTable.hpp>
#include <sgpp/optimization/test_problems/unconstrained/IncreasingPower.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Michalewicz.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Mladineo.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Perm.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rosenbrock.hpp>
#include <sgpp/optimization/test_problems/unconstrained/SHCB.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Schwefel06.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Schwefel22.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Schwefel26.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>
#include <sgpp/optimization/test_problems/unconstrained/TremblingParabola.hpp>

#include <sgpp/optimization/test_problems/constrained/Floudas.hpp>
#include <sgpp/optimization/test_problems/constrained/G03.hpp>
#include <sgpp/optimization/test_problems/constrained/G04.hpp>
#include <sgpp/optimization/test_problems/constrained/G04Squared.hpp>
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

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <vector>

using sgpp::base::Printer;
using sgpp::base::RandomNumberGenerator;
using sgpp::base::ScalarFunction;
using sgpp::base::VectorFunction;
using sgpp::optimization::test_problems::ConstrainedTestProblem;
using sgpp::optimization::test_problems::UnconstrainedTestProblem;

BOOST_AUTO_TEST_CASE(TestUnconstrainedTestProblem) {
  // Test unconstrained test problems in sgpp::optimization::test_problems.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 6;
  std::vector<std::unique_ptr<UnconstrainedTestProblem>> testProblems;
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::AbsoluteValue(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Ackley(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Alpine02(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Beale()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Branin01()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Branin02()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::BubbleWrap(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::EasomYang(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Eggholder()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::GoldsteinPrice()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Griewank(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Hartman3()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Hartman6()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Himmelblau()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::HoelderTable()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::IncreasingPower(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Michalewicz()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Mladineo()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Perm(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Rastrigin(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Rosenbrock(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Schwefel06()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Schwefel22(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Schwefel26(d)));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::SHCB()));
  testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
      new sgpp::optimization::test_problems::Sphere(d)));

  for (size_t p = 1; p <= 5; p++) {
    testProblems.push_back(std::unique_ptr<UnconstrainedTestProblem>(
        new sgpp::optimization::test_problems::TremblingParabola(d, p)));
  }

  for (const auto& problem : testProblems) {
    sgpp::optimization::test_problems::TestScalarFunction& f = problem->getObjectiveFunction();
    const size_t d = f.getNumberOfParameters();

    // test cloning
    std::unique_ptr<ScalarFunction> fClone(nullptr);
    f.clone(fClone);

    // check displacement
    sgpp::base::DataVector displacement(d, 0.42);
    f.setDisplacement(displacement);
    sgpp::base::DataVector displacement2(f.getDisplacement());

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(displacement[t], displacement2[t]);
    }

    problem->generateDisplacement();
    displacement = problem->getDisplacement();
    displacement2 = f.getDisplacement();

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(displacement[t], displacement2[t]);
    }

    displacement.setAll(0.42);
    problem->setDisplacement(displacement);
    displacement2 = problem->getDisplacement();

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(displacement[t], displacement2[t]);
    }

    // displace function randomly
    problem->generateDisplacement();

    // test minimal point
    sgpp::base::DataVector x(d), xOpt(0);
    const double fOpt = problem->getOptimalPoint(xOpt);

    BOOST_CHECK_EQUAL(xOpt.getSize(), d);

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_GE(xOpt[t], 0.0);
      BOOST_CHECK_LE(xOpt[t], 1.0);
    }

    BOOST_CHECK_SMALL(fOpt - f.eval(xOpt), 1e-12);

    // test if xOpt is minimal point for a sample of random points
    for (size_t i = 0; i < 1000; i++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = RandomNumberGenerator::getInstance().getUniformRN();
      }

      // use cloned function to test the cloning
      BOOST_CHECK_GE(fClone->eval(x), fOpt);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestConstrainedTestProblem) {
  // Test constrained test problems in sgpp::optimization::test_problems.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 6;
  std::vector<std::unique_ptr<ConstrainedTestProblem>> testProblems;
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::Floudas()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G03(d)));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G04()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G04Squared()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G05()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G06()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G08()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G09()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G10()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G11()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G12()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::G13()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::Simionescu()));
  testProblems.push_back(
      std::unique_ptr<ConstrainedTestProblem>(new sgpp::optimization::test_problems::Soland()));

  for (const auto& problem : testProblems) {
    sgpp::optimization::test_problems::TestScalarFunction& f = problem->getObjectiveFunction();
    sgpp::optimization::test_problems::TestVectorFunction& g =
        problem->getInequalityConstraintFunction();
    sgpp::optimization::test_problems::TestVectorFunction& h =
        problem->getEqualityConstraintFunction();
    const size_t d = f.getNumberOfParameters();

    // test cloning
    std::unique_ptr<ScalarFunction> fClone(nullptr);
    f.clone(fClone);
    std::unique_ptr<VectorFunction> gClone(nullptr);
    g.clone(gClone);
    std::unique_ptr<VectorFunction> hClone(nullptr);
    h.clone(hClone);

    // check displacement
    sgpp::base::DataVector displacement(d, 0.42);
    f.setDisplacement(displacement);
    g.setDisplacement(displacement);
    h.setDisplacement(displacement);
    sgpp::base::DataVector displacement2f(f.getDisplacement());
    sgpp::base::DataVector displacement2g(g.getDisplacement());
    sgpp::base::DataVector displacement2h(h.getDisplacement());

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(displacement[t], displacement2f[t]);
      BOOST_CHECK_EQUAL(displacement[t], displacement2g[t]);
      BOOST_CHECK_EQUAL(displacement[t], displacement2h[t]);
    }

    problem->generateDisplacement();
    displacement = problem->getDisplacement();
    displacement2f = f.getDisplacement();
    displacement2g = g.getDisplacement();
    displacement2h = h.getDisplacement();

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(displacement[t], displacement2f[t]);
      BOOST_CHECK_EQUAL(displacement[t], displacement2g[t]);
      BOOST_CHECK_EQUAL(displacement[t], displacement2h[t]);
    }

    displacement.setAll(0.42);
    problem->setDisplacement(displacement);
    displacement2f = problem->getDisplacement();

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(displacement[t], displacement2f[t]);
    }

    // displace function randomly
    problem->generateDisplacement();

    // test minimal point
    sgpp::base::DataVector x(d), xOpt(0);
    const double fOpt = problem->getOptimalPoint(xOpt);

    BOOST_CHECK_EQUAL(xOpt.getSize(), d);

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_GE(xOpt[t], 0.0);
      BOOST_CHECK_LE(xOpt[t], 1.0);
    }

    BOOST_CHECK_SMALL(fOpt - f.eval(xOpt), 1e-6);

    const size_t mG = g.getNumberOfComponents();
    const size_t mH = h.getNumberOfComponents();
    sgpp::base::DataVector gx(mG);
    sgpp::base::DataVector hx(mH);

    g.eval(xOpt, gx);
    h.eval(xOpt, hx);

    for (size_t t = 0; t < mG; t++) {
      BOOST_CHECK_LE(gx[t], 1e-10);
    }

    for (size_t t = 0; t < mH; t++) {
      BOOST_CHECK_SMALL(hx[t], 1e-3);
    }

    // test if xOpt is minimal point for a sample of random points
    for (size_t i = 0; i < 1000; i++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = RandomNumberGenerator::getInstance().getUniformRN();
      }

      bool isFeasible = true;

      // use cloned functions to test the cloning
      gClone->eval(x, gx);
      hClone->eval(x, hx);

      for (size_t t = 0; t < mG; t++) {
        if (gx[t] > 0.0) {
          isFeasible = false;
          break;
        }
      }

      for (size_t t = 0; t < mH; t++) {
        if (hx[t] != 0.0) {
          isFeasible = false;
          break;
        }
      }

      // we don't test if the point is infeasible
      // unfortunately, there is no general method to sample only feasible points
      if (isFeasible) {
        BOOST_CHECK_GE(fClone->eval(x), fOpt);
      }
    }
  }
}

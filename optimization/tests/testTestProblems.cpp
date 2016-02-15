// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/optimization/test_problems/unconstrained/AbsoluteValue.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Ackley.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Beale.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Branin.hpp>
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
#include <sgpp/optimization/test_problems/unconstrained/Schwefel.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>
#include <sgpp/optimization/test_problems/unconstrained/TremblingParabola.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <vector>

using SGPP::optimization::Printer;
using SGPP::optimization::RandomNumberGenerator;
using SGPP::optimization::ScalarFunction;

BOOST_AUTO_TEST_CASE(TestUnconstrainedTestProblem) {
  // Test unconstrained test problems in SGPP::optimization::test_problems.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 6;
  std::vector<std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>>
      testProblems;
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::AbsoluteValue(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Ackley(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Beale())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Branin())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::BubbleWrap(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::EasomYang(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Eggholder())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::GoldsteinPrice())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Griewank(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Hartman3())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Hartman6())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Himmelblau())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::HoelderTable())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::IncreasingPower(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Michalewicz())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Mladineo())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Perm(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Rastrigin(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Rosenbrock(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Schwefel(d))));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::SHCB())));
  testProblems.push_back(
    std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                new SGPP::optimization::test_problems::Sphere(d))));

  for (size_t p = 1; p <= 5; p++) {
    testProblems.push_back(
      std::move(std::unique_ptr<SGPP::optimization::test_problems::UnconstrainedTestProblem>(
                  new SGPP::optimization::test_problems::TremblingParabola(d, p))));
  }

  for (const auto& problem : testProblems) {
    SGPP::optimization::test_problems::TestScalarFunction& f = problem->getObjectiveFunction();
    const size_t d = f.getNumberOfParameters();

    // test cloning
    std::unique_ptr<ScalarFunction> fClone(nullptr);
    f.clone(fClone);

    // check displacement
    SGPP::base::DataVector displacement(d, 0.42);
    f.setDisplacement(displacement);
    SGPP::base::DataVector displacement2(f.getDisplacement());

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
    SGPP::base::DataVector x(d), xOpt(0);
    const SGPP::float_t fOpt = problem->getOptimalPoint(xOpt);

    BOOST_CHECK_EQUAL(xOpt.getSize(), d);

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_GE(xOpt[t], 0.0);
      BOOST_CHECK_LE(xOpt[t], 1.0);
    }

#if USE_DOUBLE_PRECISION
    BOOST_CHECK_SMALL(fOpt - f.eval(xOpt), 1e-12);
#else
    BOOST_CHECK_SMALL(fOpt - f.eval(xOpt),
                      static_cast<SGPP::float_t>(1e-3));
#endif

    // test if xopt is minimal point for a sample of random points
    for (size_t i = 0; i < 1000; i++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = RandomNumberGenerator::getInstance().getUniformRN();
      }

      // use cloned function to test the cloning
      BOOST_CHECK_GE(fClone->eval(x), fOpt);
    }
  }
}

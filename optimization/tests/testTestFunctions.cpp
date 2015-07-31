#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/test/Ackley.hpp>
#include <sgpp/optimization/function/test/Beale.hpp>
#include <sgpp/optimization/function/test/Branin.hpp>
#include <sgpp/optimization/function/test/Easom.hpp>
#include <sgpp/optimization/function/test/Eggholder.hpp>
#include <sgpp/optimization/function/test/GoldsteinPrice.hpp>
#include <sgpp/optimization/function/test/Griewank.hpp>
#include <sgpp/optimization/function/test/Hartman3.hpp>
#include <sgpp/optimization/function/test/Hartman6.hpp>
#include <sgpp/optimization/function/test/Himmelblau.hpp>
#include <sgpp/optimization/function/test/HoelderTable.hpp>
#include <sgpp/optimization/function/test/Michalewicz.hpp>
#include <sgpp/optimization/function/test/Mladineo.hpp>
#include <sgpp/optimization/function/test/Rastrigin.hpp>
#include <sgpp/optimization/function/test/Rosenbrock.hpp>
#include <sgpp/optimization/function/test/SHCB.hpp>
#include <sgpp/optimization/function/test/Schwefel.hpp>
#include <sgpp/optimization/function/test/Sphere.hpp>
#include <sgpp/optimization/function/test/TestFunction.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

using namespace SGPP;
using namespace SGPP::optimization;

BOOST_AUTO_TEST_CASE(TestFunctions) {
  // Test SGPP::optimization::test_functions::TestFunction.
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t d = 10;
  std::vector<std::unique_ptr<test_functions::TestFunction>> testFunctions;
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Ackley(d))));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Beale())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Branin())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Easom())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Eggholder())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::GoldsteinPrice())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Griewank(d))));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Hartman3())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Hartman6())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Himmelblau())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::HoelderTable())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Michalewicz())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Mladineo())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Rastrigin(d))));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Rosenbrock(d))));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Schwefel(d))));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::SHCB())));
  testFunctions.push_back(
    std::move(std::unique_ptr<test_functions::TestFunction>(
                new test_functions::Sphere(d))));

  for (const auto& fcn : testFunctions) {
    const size_t d = fcn->getDimension();
    // displace function randomly
    fcn->generateDisplacement();

    // generate random point, displace and reverse
    base::DataVector xl(d);

    for (size_t t = 0; t < d; t++) {
      xl[t] = randomNumberGenerator.getUniformRN();
    }

    base::DataVector x(xl);
    fcn->displaceVector(x);
    const SGPP::float_t f1 = fcn->evalUndisplaced(x);
    fcn->reverseDisplaceVector(x);
    const SGPP::float_t f2 = fcn->eval(x);

    // test displaceVector/reverseDisplaceVector
    for (size_t t = 0; t < d; t++) {
#if USE_DOUBLE_PRECISION == 1
      BOOST_CHECK_CLOSE(xl[t], x[t], 1e-10);
#else
      BOOST_CHECK_CLOSE(xl[t], x[t], 1e-5);
#endif
    }

    // test eval/evalUndisplaced
    BOOST_CHECK_EQUAL(f1, f2);

    base::DataVector xOpt(0);
    const SGPP::float_t fOpt = fcn->getOptimalPoint(xOpt);

    // test minimal point
    BOOST_CHECK_EQUAL(xOpt.getSize(), d);

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_GE(xOpt[t], 0.0);
      BOOST_CHECK_LE(xOpt[t], 1.0);
    }

#if USE_DOUBLE_PRECISION == 1
    BOOST_CHECK_EQUAL(fOpt, fcn->eval(xOpt) );
#else
    BOOST_CHECK_SMALL( std::abs( fOpt-fcn->eval(xOpt) ), float_t(1e-6) );
#endif

    // test if xopt is minimal point for a sample of random points
    for (size_t i = 0; i < 1000; i++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = randomNumberGenerator.getUniformRN();
      }

      BOOST_CHECK_GE(fcn->eval(x), fOpt);
    }
  }
}

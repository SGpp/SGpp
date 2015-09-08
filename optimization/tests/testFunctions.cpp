#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/scalar/ScalarComponent.hpp>
#include <sgpp/optimization/function/scalar/ScalarComponentGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarComponentHessian.hpp>

#include <sgpp/optimization/function/scalar/test/Ackley.hpp>
#include <sgpp/optimization/function/scalar/test/Beale.hpp>
#include <sgpp/optimization/function/scalar/test/Branin.hpp>
#include <sgpp/optimization/function/scalar/test/Easom.hpp>
#include <sgpp/optimization/function/scalar/test/Eggholder.hpp>
#include <sgpp/optimization/function/scalar/test/GoldsteinPrice.hpp>
#include <sgpp/optimization/function/scalar/test/Griewank.hpp>
#include <sgpp/optimization/function/scalar/test/Hartman3.hpp>
#include <sgpp/optimization/function/scalar/test/Hartman6.hpp>
#include <sgpp/optimization/function/scalar/test/Himmelblau.hpp>
#include <sgpp/optimization/function/scalar/test/HoelderTable.hpp>
#include <sgpp/optimization/function/scalar/test/Michalewicz.hpp>
#include <sgpp/optimization/function/scalar/test/Mladineo.hpp>
#include <sgpp/optimization/function/scalar/test/Rastrigin.hpp>
#include <sgpp/optimization/function/scalar/test/Rosenbrock.hpp>
#include <sgpp/optimization/function/scalar/test/SHCB.hpp>
#include <sgpp/optimization/function/scalar/test/Schwefel.hpp>
#include <sgpp/optimization/function/scalar/test/Sphere.hpp>
#include <sgpp/optimization/function/scalar/test/TestFunction.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

const bool use_double_precision =
#if USE_DOUBLE_PRECISION
  true;
#else
  false;
#endif /* USE_DOUBLE_PRECISION */

using namespace SGPP;
using namespace SGPP::optimization;

class ScalarTestFunction : public ScalarFunction {
  public:
    ScalarTestFunction(size_t d) : ScalarFunction(d) {}

    SGPP::float_t eval(const base::DataVector& x) {
      return x.sum();
    }

    virtual void clone(std::unique_ptr<ScalarFunction>& clone) const {
      clone = std::unique_ptr<ScalarFunction>(
                new ScalarTestFunction(*this));
    }
};

class ScalarTestGradient : public ScalarFunctionGradient {
  public:
    ScalarTestGradient(size_t d) : ScalarFunctionGradient(d) {}

    SGPP::float_t eval(const base::DataVector& x,
                       base::DataVector& gradient) {
      for (size_t t = 0; t < d; t++) {
        gradient[t] = static_cast<SGPP::float_t>(t) * x[t];
      }

      return x.sum();
    }

    virtual void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const {
      clone = std::unique_ptr<ScalarFunctionGradient>(
                new ScalarTestGradient(*this));
    }
};

class ScalarTestHessian : public ScalarFunctionHessian {
  public:
    ScalarTestHessian(size_t d) : ScalarFunctionHessian(d) {}

    SGPP::float_t eval(const base::DataVector& x,
                       base::DataVector& gradient,
                       base::DataMatrix& hessian) {
      for (size_t t = 0; t < d; t++) {
        gradient[t] = static_cast<SGPP::float_t>(t) * x[t];

        for (size_t t2 = 0; t2 < d; t2++) {
          hessian(t, t2) = static_cast<SGPP::float_t>(t) * x[t] +
                           static_cast<SGPP::float_t>(t2) * x[t2];
        }
      }

      return x.sum();
    }

    virtual void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const {
      clone = std::unique_ptr<ScalarFunctionHessian>(
                new ScalarTestHessian(*this));
    }
};

class VectorTestFunction : public VectorFunction {
  public:
    VectorTestFunction(size_t d, size_t m) : VectorFunction(d, m) {}

    void eval(const base::DataVector& x,
              base::DataVector& value) {
      for (size_t i = 0; i < m; i++) {
        value[i] = static_cast<SGPP::float_t>(i) * x.sum();
      }
    }

    virtual void clone(std::unique_ptr<VectorFunction>& clone) const {
      clone = std::unique_ptr<VectorFunction>(
                new VectorTestFunction(*this));
    }
};

class VectorTestGradient : public VectorFunctionGradient {
  public:
    VectorTestGradient(size_t d, size_t m) : VectorFunctionGradient(d, m) {}

    void eval(const base::DataVector& x,
              base::DataVector& value,
              base::DataMatrix& gradient) {
      for (size_t i = 0; i < m; i++) {
        value[i] = static_cast<SGPP::float_t>(i) * x.sum();

        for (size_t t = 0; t < d; t++) {
          gradient(i, t) = static_cast<SGPP::float_t>(i) *
                           static_cast<SGPP::float_t>(t) * x[t];
        }
      }
    }

    virtual void clone(std::unique_ptr<VectorFunctionGradient>& clone) const {
      clone = std::unique_ptr<VectorFunctionGradient>(
                new VectorTestGradient(*this));
    }
};

class VectorTestHessian : public VectorFunctionHessian {
  public:
    VectorTestHessian(size_t d, size_t m) : VectorFunctionHessian(d, m) {}

    void eval(const base::DataVector& x,
              base::DataVector& value,
              base::DataMatrix& gradient,
              std::vector<base::DataMatrix>& hessian) {
      for (size_t i = 0; i < m; i++) {
        value[i] = static_cast<SGPP::float_t>(i) * x.sum();

        for (size_t t = 0; t < d; t++) {
          gradient(i, t) = static_cast<SGPP::float_t>(i) *
                           static_cast<SGPP::float_t>(t) * x[t];

          for (size_t t2 = 0; t2 < d; t2++) {
            hessian[i](t, t2) = static_cast<SGPP::float_t>(i) *
                                static_cast<SGPP::float_t>(t) * x[t] *
                                static_cast<SGPP::float_t>(t2) * x[t];
          }
        }
      }
    }

    virtual void clone(std::unique_ptr<VectorFunctionHessian>& clone) const {
      clone = std::unique_ptr<VectorFunctionHessian>(
                new VectorTestHessian(*this));
    }
};

BOOST_AUTO_TEST_CASE(TestScalarComponent) {
  // Test SGPP::optimization::ScalarComponent for scalar functions.
  base::DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  base::DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestFunction f(3);
    ScalarComponent g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunction> g2;
    g.clone(g2);
    BOOST_CHECK_EQUAL(f.eval(x), g2->eval(y));

    BOOST_CHECK_THROW(ScalarComponent(f, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ScalarComponent(f, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }

  {
    VectorTestFunction f(3, 4);
    base::DataVector fx(4);
    f.eval(x, fx);
    ScalarComponent g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunction> g2;
    g.clone(g2);
    BOOST_CHECK_EQUAL(fx[1], g2->eval(y));

    BOOST_CHECK_THROW(ScalarComponent(f, 1, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ScalarComponent(f, 1, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(TestScalarComponentGradient) {
  // Test SGPP::optimization::ScalarComponentGradient for scalar functions.
  base::DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  base::DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestGradient f(3);
    ScalarComponentGradient g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionGradient> g2;
    g.clone(g2);
    base::DataVector gradF(3), gradG(2);
    BOOST_CHECK_EQUAL(f.eval(x, gradF), g2->eval(y, gradG));
    BOOST_CHECK_EQUAL(gradF[0], gradG[0]);
    BOOST_CHECK_EQUAL(gradF[2], gradG[1]);

    BOOST_CHECK_THROW(ScalarComponentGradient(f, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ScalarComponentGradient(f, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }

  {
    VectorTestGradient f(3, 4);
    base::DataVector fx(4);
    base::DataMatrix gradF(4, 3);
    f.eval(x, fx, gradF);
    ScalarComponentGradient g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionGradient> g2;
    g.clone(g2);
    base::DataVector gradG(2);
    BOOST_CHECK_EQUAL(fx[1], g2->eval(y, gradG));
    BOOST_CHECK_EQUAL(gradF(1, 0), gradG[0]);
    BOOST_CHECK_EQUAL(gradF(1, 2), gradG[1]);

    BOOST_CHECK_THROW(ScalarComponentGradient(f, 1, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ScalarComponentGradient(f, 1, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(TestScalarComponentHessian) {
  // Test SGPP::optimization::ScalarComponentHessian for scalar functions.
  base::DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  base::DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestHessian f(3);
    ScalarComponentHessian g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionHessian> g2;
    g.clone(g2);
    base::DataVector gradF(3), gradG(2);
    base::DataMatrix hessF(3, 3), hessG(2, 2);
    BOOST_CHECK_EQUAL(f.eval(x, gradF, hessF), g2->eval(y, gradG, hessG));
    BOOST_CHECK_EQUAL(gradF[0], gradG[0]);
    BOOST_CHECK_EQUAL(gradF[2], gradG[1]);
    BOOST_CHECK_EQUAL(hessF(0, 0), hessG(0, 0));
    BOOST_CHECK_EQUAL(hessF(0, 2), hessG(0, 1));
    BOOST_CHECK_EQUAL(hessF(2, 0), hessG(1, 0));
    BOOST_CHECK_EQUAL(hessF(2, 2), hessG(1, 1));

    BOOST_CHECK_THROW(ScalarComponentHessian(f, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ScalarComponentHessian(f, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }

  {
    VectorTestHessian f(3, 4);
    base::DataVector fx(4);
    base::DataMatrix gradF(4, 3);
    std::vector<base::DataMatrix> hessF(4, base::DataMatrix(3, 3));
    f.eval(x, fx, gradF, hessF);
    ScalarComponentHessian g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionHessian> g2;
    g.clone(g2);
    base::DataVector gradG(2);
    base::DataMatrix hessG(2, 2);
    BOOST_CHECK_EQUAL(fx[1], g2->eval(y, gradG, hessG));
    BOOST_CHECK_EQUAL(gradF(1, 0), gradG[0]);
    BOOST_CHECK_EQUAL(gradF(1, 2), gradG[1]);
    BOOST_CHECK_EQUAL(hessF[1](0, 0), hessG(0, 0));
    BOOST_CHECK_EQUAL(hessF[1](0, 2), hessG(0, 1));
    BOOST_CHECK_EQUAL(hessF[1](2, 0), hessG(1, 0));
    BOOST_CHECK_EQUAL(hessF[1](2, 2), hessG(1, 1));

    BOOST_CHECK_THROW(ScalarComponentHessian(f, 1, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ScalarComponentHessian(f, 1, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(TestTestFunctions) {
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

    // test cloning
    std::unique_ptr<ScalarFunction> fcn2(nullptr);
    fcn->clone(fcn2);

    // check displacement
    const SGPP::float_t stdDev = 0.01;
    fcn->generateDisplacement(stdDev);
    BOOST_CHECK_EQUAL(fcn->getStandardDeviation(), stdDev);
    base::DataVector displacement(0);
    fcn->getDisplacement(displacement);

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
      BOOST_CHECK_CLOSE(xl[t], x[t], (use_double_precision ? 1e-10 : 1e-5));
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
    BOOST_CHECK_EQUAL(fOpt, fcn->eval(xOpt));
#else
    BOOST_CHECK_SMALL(fOpt - fcn->eval(xOpt),
                      static_cast<SGPP::float_t>(1e-6));
#endif

    // test if xopt is minimal point for a sample of random points
    for (size_t i = 0; i < 1000; i++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = randomNumberGenerator.getUniformRN();
      }

      // use cloned function to test the cloning
      BOOST_CHECK_GE(fcn2->eval(x), fOpt);
    }
  }
}

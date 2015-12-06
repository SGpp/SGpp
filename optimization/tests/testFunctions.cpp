#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/scalar/ComponentScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ComponentScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ComponentScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionHessian.hpp>

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

BOOST_AUTO_TEST_CASE(TestComponentScalarFunction) {
  // Test SGPP::optimization::ComponentScalarFunction.
  base::DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  base::DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestFunction f(3);
    ComponentScalarFunction g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunction> g2;
    g.clone(g2);
    BOOST_CHECK_EQUAL(f.eval(x), g2->eval(y));

    BOOST_CHECK_THROW(ComponentScalarFunction(f, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ComponentScalarFunction(f, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }

  {
    VectorTestFunction f(3, 4);
    base::DataVector fx(4);
    f.eval(x, fx);
    ComponentScalarFunction g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunction> g2;
    g.clone(g2);
    BOOST_CHECK_EQUAL(fx[1], g2->eval(y));

    BOOST_CHECK_THROW(ComponentScalarFunction(f, 1, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ComponentScalarFunction(f, 1, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(TestComponentScalarFunctionGradient) {
  // Test SGPP::optimization::ComponentScalarFunctionGradient.
  base::DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  base::DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestGradient f(3);
    ComponentScalarFunctionGradient g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionGradient> g2;
    g.clone(g2);
    base::DataVector gradF(3), gradG(2);
    BOOST_CHECK_EQUAL(f.eval(x, gradF), g2->eval(y, gradG));
    BOOST_CHECK_EQUAL(gradF[0], gradG[0]);
    BOOST_CHECK_EQUAL(gradF[2], gradG[1]);

    BOOST_CHECK_THROW(ComponentScalarFunctionGradient(f, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ComponentScalarFunctionGradient(f, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }

  {
    VectorTestGradient f(3, 4);
    base::DataVector fx(4);
    base::DataMatrix gradF(4, 3);
    f.eval(x, fx, gradF);
    ComponentScalarFunctionGradient g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionGradient> g2;
    g.clone(g2);
    base::DataVector gradG(2);
    BOOST_CHECK_EQUAL(fx[1], g2->eval(y, gradG));
    BOOST_CHECK_EQUAL(gradF(1, 0), gradG[0]);
    BOOST_CHECK_EQUAL(gradF(1, 2), gradG[1]);

    BOOST_CHECK_THROW(ComponentScalarFunctionGradient(f, 1, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ComponentScalarFunctionGradient(f, 1, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(TestComponentScalarFunctionHessian) {
  // Test SGPP::optimization::ComponentScalarFunctionHessian.
  base::DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  base::DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestHessian f(3);
    ComponentScalarFunctionHessian g(f, {NAN, 0.34, NAN});
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

    BOOST_CHECK_THROW(ComponentScalarFunctionHessian(f, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ComponentScalarFunctionHessian(f, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }

  {
    VectorTestHessian f(3, 4);
    base::DataVector fx(4);
    base::DataMatrix gradF(4, 3);
    std::vector<base::DataMatrix> hessF(4, base::DataMatrix(3, 3));
    f.eval(x, fx, gradF, hessF);
    ComponentScalarFunctionHessian g(f, 1, {NAN, 0.34, NAN});
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

    BOOST_CHECK_THROW(ComponentScalarFunctionHessian(f, 1, {NAN, NAN}),
                      std::runtime_error);
    BOOST_CHECK_THROW(ComponentScalarFunctionHessian(f, 1, {NAN, NAN, NAN, NAN}),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(TestWrapperScalarFunction) {
  // Test SGPP::optimization::TestWrapperScalarFunction.
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 3;
  const size_t N = 100;
  base::DataVector x(d);
  ScalarTestFunction f1(d);
  WrapperScalarFunction f2(d, [](const base::DataVector & x) {
    return x.sum();
  });
  std::unique_ptr<ScalarFunction> f2Clone;
  f2.clone(f2Clone);

  for (size_t i = 0; i < N; i++) {
    for (size_t t = 0; t < d; t++) {
      x[t] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    BOOST_CHECK_EQUAL(f1.eval(x), f2Clone->eval(x));
  }
}

BOOST_AUTO_TEST_CASE(TestWrapperScalarFunctionGradient) {
  // Test SGPP::optimization::TestWrapperScalarFunctionGradient.
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 3;
  const size_t N = 100;
  base::DataVector x(d);
  base::DataVector gradient1(d), gradient2(d);
  ScalarTestGradient f1(d);
  WrapperScalarFunctionGradient f2(d, [](const base::DataVector & x,
  base::DataVector & gradient) {
    for (size_t t = 0; t < d; t++) {
      gradient[t] = static_cast<SGPP::float_t>(t) * x[t];
    }

    return x.sum();
  });
  std::unique_ptr<ScalarFunctionGradient> f2Clone;
  f2.clone(f2Clone);

  for (size_t i = 0; i < N; i++) {
    for (size_t t = 0; t < d; t++) {
      x[t] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    BOOST_CHECK_EQUAL(f1.eval(x, gradient1), f2Clone->eval(x, gradient2));
    BOOST_CHECK_EQUAL(gradient2.getSize(), d);

    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(gradient1[t], gradient2[t]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestWrapperScalarFunctionHessian) {
  // Test SGPP::optimization::TestWrapperScalarFunctionHessian.
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 3;
  const size_t N = 100;
  base::DataVector x(d);
  base::DataVector gradient1(d), gradient2(d);
  base::DataMatrix hessian1(d, d), hessian2(d, d);
  ScalarTestHessian f1(d);
  WrapperScalarFunctionHessian f2(d, [](const base::DataVector & x,
                                        base::DataVector & gradient,
  base::DataMatrix & hessian) {
    for (size_t t = 0; t < d; t++) {
      gradient[t] = static_cast<SGPP::float_t>(t) * x[t];

      for (size_t t2 = 0; t2 < d; t2++) {
        hessian(t, t2) = static_cast<SGPP::float_t>(t) * x[t] +
                         static_cast<SGPP::float_t>(t2) * x[t2];
      }
    }

    return x.sum();
  });
  std::unique_ptr<ScalarFunctionHessian> f2Clone;
  f2.clone(f2Clone);

  for (size_t i = 0; i < N; i++) {
    for (size_t t = 0; t < d; t++) {
      x[t] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    BOOST_CHECK_EQUAL(f1.eval(x, gradient1, hessian1),
                      f2Clone->eval(x, gradient2, hessian2));
    BOOST_CHECK_EQUAL(gradient2.getSize(), d);
    BOOST_CHECK_EQUAL(hessian2.getNrows(), d);
    BOOST_CHECK_EQUAL(hessian2.getNcols(), d);

    for (size_t t1 = 0; t1 < d; t1++) {
      for (size_t t2 = 0; t2 < d; t2++) {
        BOOST_CHECK_EQUAL(hessian1(t1, t2), hessian2(t1, t2));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestWrapperVectorFunction) {
  // Test SGPP::optimization::TestWrapperVectorFunction.
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 3;
  const size_t m = 4;
  const size_t N = 100;
  base::DataVector x(d);
  base::DataVector value1(m), value2(m);
  VectorTestFunction f1(d, m);
  WrapperVectorFunction f2(d, m, [](const base::DataVector & x,
  base::DataVector & value) {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<SGPP::float_t>(i) * x.sum();
    }
  });
  std::unique_ptr<VectorFunction> f2Clone;
  f2.clone(f2Clone);

  for (size_t i = 0; i < N; i++) {
    for (size_t t = 0; t < d; t++) {
      x[t] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    f1.eval(x, value1);
    f2Clone->eval(x, value2);

    BOOST_CHECK_EQUAL(value2.getSize(), m);

    for (size_t j = 0; j < m; j++) {
      BOOST_CHECK_EQUAL(value1[j], value2[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestWrapperVectorFunctionGradient) {
  // Test SGPP::optimization::TestWrapperVectorFunctionGradient.
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 3;
  const size_t m = 4;
  const size_t N = 100;
  base::DataVector x(d);
  base::DataVector value1(m), value2(m);
  base::DataMatrix gradient1(m, d), gradient2(m, d);
  VectorTestGradient f1(d, m);
  WrapperVectorFunctionGradient f2(d, m, [](const base::DataVector & x,
                                   base::DataVector & value,
  base::DataMatrix & gradient) {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<SGPP::float_t>(i) * x.sum();

      for (size_t t = 0; t < d; t++) {
        gradient(i, t) = static_cast<SGPP::float_t>(i) *
                         static_cast<SGPP::float_t>(t) * x[t];
      }
    }
  });
  std::unique_ptr<VectorFunctionGradient> f2Clone;
  f2.clone(f2Clone);

  for (size_t i = 0; i < N; i++) {
    for (size_t t = 0; t < d; t++) {
      x[t] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    f1.eval(x, value1, gradient1);
    f2Clone->eval(x, value2, gradient2);

    BOOST_CHECK_EQUAL(value2.getSize(), m);
    BOOST_CHECK_EQUAL(gradient2.getNrows(), m);
    BOOST_CHECK_EQUAL(gradient2.getNcols(), d);

    for (size_t j = 0; j < m; j++) {
      BOOST_CHECK_EQUAL(value1[j], value2[j]);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_EQUAL(gradient1(j, t), gradient2(j, t));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestWrapperVectorFunctionHessian) {
  // Test SGPP::optimization::TestWrapperVectorFunctionHessian.
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 3;
  const size_t m = 4;
  const size_t N = 100;
  base::DataVector x(d);
  base::DataVector value1(m), value2(m);
  base::DataMatrix gradient1(m, d), gradient2(m, d);
  std::vector<base::DataMatrix> hessian1(m, base::DataMatrix(d, d)),
      hessian2(m, base::DataMatrix(d, d));
  VectorTestHessian f1(d, m);
  WrapperVectorFunctionHessian f2(d, m, [](const base::DataVector & x,
                                  base::DataVector & value,
                                  base::DataMatrix & gradient,
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
  });
  std::unique_ptr<VectorFunctionHessian> f2Clone;
  f2.clone(f2Clone);

  for (size_t i = 0; i < N; i++) {
    for (size_t t = 0; t < d; t++) {
      x[t] = RandomNumberGenerator::getInstance().getUniformRN();
    }


    f1.eval(x, value1, gradient1, hessian1);
    f2Clone->eval(x, value2, gradient2, hessian2);

    BOOST_CHECK_EQUAL(value2.getSize(), m);
    BOOST_CHECK_EQUAL(gradient2.getNrows(), m);
    BOOST_CHECK_EQUAL(gradient2.getNcols(), d);
    BOOST_CHECK_EQUAL(hessian2.size(), m);

    for (size_t j = 0; j < m; j++) {
      BOOST_CHECK_EQUAL(value1[j], value2[j]);
      BOOST_CHECK_EQUAL(hessian2[j].getNrows(), d);
      BOOST_CHECK_EQUAL(hessian2[j].getNcols(), d);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_EQUAL(gradient1(j, t), gradient2(j, t));

        for (size_t t2 = 0; t2 < d; t2++) {
          BOOST_CHECK_EQUAL(hessian1[j](t, t2), hessian2[j](t, t2));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestTestFunctions) {
  // Test SGPP::optimization::test_functions::TestFunction.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

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
    const size_t d = fcn->getNumberOfParameters();

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
      xl[t] = RandomNumberGenerator::getInstance().getUniformRN();
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

#if USE_DOUBLE_PRECISION
    BOOST_CHECK_SMALL(fOpt - fcn->eval(xOpt), 1e-10);
#else
    BOOST_CHECK_SMALL(fOpt - fcn->eval(xOpt), 1e-6f);
#endif /* USE_DOUBLE_PRECISION */

    // test if xopt is minimal point for a sample of random points
    for (size_t i = 0; i < 1000; i++) {
      for (size_t t = 0; t < d; t++) {
        x[t] = RandomNumberGenerator::getInstance().getUniformRN();
      }

      // use cloned function to test the cloning
      BOOST_CHECK_GE(fcn2->eval(x), fOpt);
    }
  }
}

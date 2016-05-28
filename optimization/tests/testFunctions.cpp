// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/optimization/function/scalar/ComponentScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ComponentScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ComponentScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionHessian.hpp>

#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include <vector>

#include "CheckEqualFunction.hpp"

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::optimization::ComponentScalarFunction;
using sgpp::optimization::ComponentScalarFunctionGradient;
using sgpp::optimization::ComponentScalarFunctionHessian;
using sgpp::optimization::RandomNumberGenerator;
using sgpp::optimization::ScalarFunction;
using sgpp::optimization::ScalarFunctionGradient;
using sgpp::optimization::ScalarFunctionHessian;
using sgpp::optimization::VectorFunction;
using sgpp::optimization::VectorFunctionGradient;
using sgpp::optimization::VectorFunctionHessian;
using sgpp::optimization::WrapperScalarFunction;
using sgpp::optimization::WrapperScalarFunctionGradient;
using sgpp::optimization::WrapperScalarFunctionHessian;
using sgpp::optimization::WrapperVectorFunction;
using sgpp::optimization::WrapperVectorFunctionGradient;
using sgpp::optimization::WrapperVectorFunctionHessian;

class ScalarTestFunction : public ScalarFunction {
 public:
  explicit ScalarTestFunction(size_t d) : ScalarFunction(d) {}

  ~ScalarTestFunction() override {}

  double eval(const DataVector& x) override {
    return x.sum();
  }

  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(
              new ScalarTestFunction(*this));
  }
};

class ScalarTestGradient : public ScalarFunctionGradient {
 public:
  explicit ScalarTestGradient(size_t d) : ScalarFunctionGradient(d) {}

  ~ScalarTestGradient() override {}

  double eval(const DataVector& x,
                     DataVector& gradient) override {
    for (size_t t = 0; t < d; t++) {
      gradient[t] = static_cast<double>(t) * x[t];
    }

    return x.sum();
  }

  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const
  override {
    clone = std::unique_ptr<ScalarFunctionGradient>(
              new ScalarTestGradient(*this));
  }
};

class ScalarTestHessian : public ScalarFunctionHessian {
 public:
  explicit ScalarTestHessian(size_t d) : ScalarFunctionHessian(d) {}

  ~ScalarTestHessian() override {}

  double eval(const DataVector& x,
                     DataVector& gradient,
                     DataMatrix& hessian) override {
    for (size_t t = 0; t < d; t++) {
      gradient[t] = static_cast<double>(t) * x[t];

      for (size_t t2 = 0; t2 < d; t2++) {
        hessian(t, t2) = static_cast<double>(t) * x[t] +
                         static_cast<double>(t2) * x[t2];
      }
    }

    return x.sum();
  }

  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const
  override {
    clone = std::unique_ptr<ScalarFunctionHessian>(
              new ScalarTestHessian(*this));
  }
};

class VectorTestFunction : public VectorFunction {
 public:
  VectorTestFunction(size_t d, size_t m) : VectorFunction(d, m) {}

  ~VectorTestFunction() override {}

  void eval(const DataVector& x,
            DataVector& value) override {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<double>(i) * x.sum();
    }
  }

  void clone(std::unique_ptr<VectorFunction>& clone) const override {
    clone = std::unique_ptr<VectorFunction>(
              new VectorTestFunction(*this));
  }
};

class VectorTestGradient : public VectorFunctionGradient {
 public:
  VectorTestGradient(size_t d, size_t m) : VectorFunctionGradient(d, m) {}

  ~VectorTestGradient() override {}

  void eval(const DataVector& x,
            DataVector& value,
            DataMatrix& gradient) override {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<double>(i) * x.sum();

      for (size_t t = 0; t < d; t++) {
        gradient(i, t) = static_cast<double>(i) *
                         static_cast<double>(t) * x[t];
      }
    }
  }

  void clone(std::unique_ptr<VectorFunctionGradient>& clone) const
  override {
    clone = std::unique_ptr<VectorFunctionGradient>(
              new VectorTestGradient(*this));
  }
};

class VectorTestHessian : public VectorFunctionHessian {
 public:
  VectorTestHessian(size_t d, size_t m) : VectorFunctionHessian(d, m) {}

  ~VectorTestHessian() override {}

  void eval(const DataVector& x,
            DataVector& value,
            DataMatrix& gradient,
            std::vector<DataMatrix>& hessian) override {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<double>(i) * x.sum();

      for (size_t t = 0; t < d; t++) {
        gradient(i, t) = static_cast<double>(i) *
                         static_cast<double>(t) * x[t];

        for (size_t t2 = 0; t2 < d; t2++) {
          hessian[i](t, t2) = static_cast<double>(i) *
                              static_cast<double>(t) * x[t] *
                              static_cast<double>(t2) * x[t];
        }
      }
    }
  }

  void clone(std::unique_ptr<VectorFunctionHessian>& clone) const
  override {
    clone = std::unique_ptr<VectorFunctionHessian>(
              new VectorTestHessian(*this));
  }
};

BOOST_AUTO_TEST_CASE(TestComponentScalarFunction) {
  // Test sgpp::optimization::ComponentScalarFunction.
  DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  DataVector y(2);
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
    DataVector fx(4);
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
  // Test sgpp::optimization::ComponentScalarFunctionGradient.
  DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestGradient f(3);
    ComponentScalarFunctionGradient g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionGradient> g2;
    g.clone(g2);
    DataVector gradF(3), gradG(2);
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
    DataVector fx(4);
    DataMatrix gradF(4, 3);
    f.eval(x, fx, gradF);
    ComponentScalarFunctionGradient g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionGradient> g2;
    g.clone(g2);
    DataVector gradG(2);
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
  // Test sgpp::optimization::ComponentScalarFunctionHessian.
  DataVector x(3);
  x[0] = 0.12;
  x[1] = 0.34;
  x[2] = 0.56;

  DataVector y(2);
  y[0] = 0.12;
  y[1] = 0.56;

  {
    ScalarTestHessian f(3);
    ComponentScalarFunctionHessian g(f, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionHessian> g2;
    g.clone(g2);
    DataVector gradF(3), gradG(2);
    DataMatrix hessF(3, 3), hessG(2, 2);
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
    DataVector fx(4);
    DataMatrix gradF(4, 3);
    std::vector<DataMatrix> hessF(4, DataMatrix(3, 3));
    f.eval(x, fx, gradF, hessF);
    ComponentScalarFunctionHessian g(f, 1, {NAN, 0.34, NAN});
    std::unique_ptr<ScalarFunctionHessian> g2;
    g.clone(g2);
    DataVector gradG(2);
    DataMatrix hessG(2, 2);
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
  // Test sgpp::optimization::TestWrapperScalarFunction.
  const size_t d = 3;
  ScalarTestFunction f1(d);
  WrapperScalarFunction f2(d, [](const DataVector & x) {
    return x.sum();
  });
  std::unique_ptr<ScalarFunction> f2Clone;
  f2.clone(f2Clone);
  checkEqualFunction(f1, *f2Clone);
}

BOOST_AUTO_TEST_CASE(TestWrapperScalarFunctionGradient) {
  // Test sgpp::optimization::TestWrapperScalarFunctionGradient.
  const size_t d = 3;
  ScalarTestGradient f1(d);
  WrapperScalarFunctionGradient f2(d, [](const DataVector & x,
  DataVector & gradient) {
    for (size_t t = 0; t < d; t++) {
      gradient[t] = static_cast<double>(t) * x[t];
    }

    return x.sum();
  });
  std::unique_ptr<ScalarFunctionGradient> f2Clone;
  f2.clone(f2Clone);
  checkEqualFunction(f1, *f2Clone);
}

BOOST_AUTO_TEST_CASE(TestWrapperScalarFunctionHessian) {
  // Test sgpp::optimization::TestWrapperScalarFunctionHessian.
  const size_t d = 3;
  ScalarTestHessian f1(d);
  WrapperScalarFunctionHessian f2(d, [](const DataVector & x,
                                        DataVector & gradient,
  DataMatrix & hessian) {
    for (size_t t = 0; t < d; t++) {
      gradient[t] = static_cast<double>(t) * x[t];

      for (size_t t2 = 0; t2 < d; t2++) {
        hessian(t, t2) = static_cast<double>(t) * x[t] +
                         static_cast<double>(t2) * x[t2];
      }
    }

    return x.sum();
  });
  std::unique_ptr<ScalarFunctionHessian> f2Clone;
  f2.clone(f2Clone);
  checkEqualFunction(f1, *f2Clone);
}

BOOST_AUTO_TEST_CASE(TestWrapperVectorFunction) {
  // Test sgpp::optimization::TestWrapperVectorFunction.
  const size_t d = 3;
  const size_t m = 4;
  DataVector value1(m), value2(m);
  VectorTestFunction f1(d, m);
  WrapperVectorFunction f2(d, m, [](const DataVector & x,
  DataVector & value) {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<double>(i) * x.sum();
    }
  });
  std::unique_ptr<VectorFunction> f2Clone;
  f2.clone(f2Clone);
  checkEqualFunction(f1, *f2Clone);
}

BOOST_AUTO_TEST_CASE(TestWrapperVectorFunctionGradient) {
  // Test sgpp::optimization::TestWrapperVectorFunctionGradient.
  const size_t d = 3;
  const size_t m = 4;
  VectorTestGradient f1(d, m);
  WrapperVectorFunctionGradient f2(d, m, [](const DataVector & x,
                                   DataVector & value,
  DataMatrix & gradient) {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<double>(i) * x.sum();

      for (size_t t = 0; t < d; t++) {
        gradient(i, t) = static_cast<double>(i) *
                         static_cast<double>(t) * x[t];
      }
    }
  });
  std::unique_ptr<VectorFunctionGradient> f2Clone;
  f2.clone(f2Clone);
  checkEqualFunction(f1, *f2Clone);
}

BOOST_AUTO_TEST_CASE(TestWrapperVectorFunctionHessian) {
  // Test sgpp::optimization::TestWrapperVectorFunctionHessian.
  const size_t d = 3;
  const size_t m = 4;
  VectorTestHessian f1(d, m);
  WrapperVectorFunctionHessian f2(d, m, [](const DataVector & x,
                                  DataVector & value,
                                  DataMatrix & gradient,
  std::vector<DataMatrix>& hessian) {
    for (size_t i = 0; i < m; i++) {
      value[i] = static_cast<double>(i) * x.sum();

      for (size_t t = 0; t < d; t++) {
        gradient(i, t) = static_cast<double>(i) *
                         static_cast<double>(t) * x[t];

        for (size_t t2 = 0; t2 < d; t2++) {
          hessian[i](t, t2) = static_cast<double>(i) *
                              static_cast<double>(t) * x[t] *
                              static_cast<double>(t2) * x[t];
        }
      }
    }
  });
  std::unique_ptr<VectorFunctionHessian> f2Clone;
  f2.clone(f2Clone);
  checkEqualFunction(f1, *f2Clone);
}

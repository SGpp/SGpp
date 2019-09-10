// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/function/vector/SparseGridResponseSurfaceBsplineVector.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::VectorFunction;

// res[t] = (1*x0 + 2*x1 + ... + d*x(d-1))^(t+1)
class multivariateTestFunction : public VectorFunction {
 public:
  multivariateTestFunction(size_t d, size_t m) : VectorFunction(d, m) {}
  ~multivariateTestFunction() override {}

  void eval(const DataVector& x, DataVector& value) override {
    for (size_t t = 0; t < m; t++) {
      double sum = 0.0;
      for (size_t s = 0; s < d; s++) {
        sum += (static_cast<double>(s) + 1) * x[s];
      }
      value[t] = std::pow(sum, t + 1);
    }
  }

  /**
   * @param x         evluation point
   * @param gradient  Jacobian of the linear combination (each row is a gradient vector)
   */
  void evalGradient(const DataVector& x, DataMatrix& gradients) {
    double sum = 0.0;
    for (size_t s = 0; s < d; s++) {
      sum += (static_cast<double>(s) + 1) * x[s];
    }
    gradients.resizeZero(m, d);
    for (size_t t = 0; t < m; t++) {
      for (size_t s = 0; s < d; s++) {
        double val = static_cast<double>(t + 1) * std::pow(sum, t) * (static_cast<double>(s) + 1);
        gradients.set(t, s, val);
      }
    }
  }

  void clone(std::unique_ptr<VectorFunction>& clone) const override {
    clone = std::unique_ptr<VectorFunction>(new multivariateTestFunction(*this));
  }

  DataVector getLowerBounds() { return DataVector(d, -2.0); }

  DataVector getUpperBounds() { return DataVector(d, 2.0); }
};

BOOST_AUTO_TEST_SUITE(testResponseSurfaces)

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorEval) {
  std::vector<double> epsilons{0.8, 1e-14, 1e-14};
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  std::vector<size_t> degrees{1, 3, 5};
  for (size_t t = 0; t < degrees.size(); t++) {
    size_t degree = degrees[t];
    DataVector lb = testFunction->getLowerBounds();
    DataVector ub = testFunction->getUpperBounds();
    sgpp::optimization::SparseGridResponseSurfaceBsplineVector reSurf(testFunction, lb, ub,
                                                                      gridType, degree);
    reSurf.regular(level);

    DataVector point(dim, 0.337);
    DataVector reSurfEval = reSurf.eval(point);
    DataVector trueEval(m);
    testFunction->eval(point, trueEval);
    trueEval.sub(reSurfEval);
    double evalErr = trueEval.sum();
    // std::cout << "evalErr = " << evalErr << "\n";
    BOOST_CHECK_SMALL(fabs(evalErr), epsilons[t]);
  }
}

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorEvalGradient) {
  double epsilon = 1e-13;
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  std::vector<size_t> degrees{3};
  for (auto& degree : degrees) {
    DataVector lb = testFunction->getLowerBounds();
    DataVector ub = testFunction->getUpperBounds();
    sgpp::optimization::SparseGridResponseSurfaceBsplineVector reSurf(testFunction, lb, ub,
                                                                      gridType, degree);
    reSurf.regular(level);

    DataVector point(dim, 1.0 / 6.0);
    DataMatrix gradients(dim, m);
    DataVector reSurfEval = reSurf.evalGradient(point, gradients);
    DataMatrix trueGradients(m, dim);
    testFunction->evalGradient(point, trueGradients);

    // std::cout << trueGradients.sum() << "\n";
    trueGradients.sub(gradients);
    double elementwiseError = trueGradients.sum();
    // std::cout << elementwiseError << "\n";
    BOOST_CHECK_SMALL(elementwiseError, epsilon);
  }
}

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorL2) {
  double epsilon = 1e-11;
  size_t dim = 3;
  size_t m = 3;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  std::vector<size_t> degrees{3};
  for (auto& degree : degrees) {
    DataVector lb = testFunction->getLowerBounds();
    DataVector ub = testFunction->getUpperBounds();
    sgpp::optimization::SparseGridResponseSurfaceBsplineVector reSurf(testFunction, lb, ub,
                                                                      gridType, degree);
    reSurf.regular(level);

    DataVector componentwiseErrors(m);
    size_t numMCPoints = 1000;
    double l2Error = reSurf.l2Error(testFunction, componentwiseErrors, numMCPoints);
    // std::cout << "l2 error: " << l2Error << "\n";
    // std::cout << "componentwise errors:\n" << componentwiseErrors.toString() << "\n";
    BOOST_CHECK_SMALL(l2Error, epsilon);
  }
}

BOOST_AUTO_TEST_SUITE_END()

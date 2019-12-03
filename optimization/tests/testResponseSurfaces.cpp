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
   * @param jacobian  Jacobian of the linear combination (each row is a gradient vector)
   */
  void evalJacobian(const DataVector& x, DataMatrix& jacobian) {
    double sum = 0.0;
    for (size_t s = 0; s < d; s++) {
      sum += (static_cast<double>(s) + 1) * x[s];
    }
    jacobian.resizeZero(m, d);
    for (size_t t = 0; t < m; t++) {
      for (size_t s = 0; s < d; s++) {
        double val = static_cast<double>(t + 1) * std::pow(sum, t) * (static_cast<double>(s) + 1);
        jacobian.set(t, s, val);
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

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorEvalJacobian) {
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
    DataMatrix jacobian(dim, m);
    DataVector reSurfEval = reSurf.evalJacobian(point, jacobian);
    DataMatrix trueJacobian(m, dim);
    testFunction->evalJacobian(point, trueJacobian);

    trueJacobian.sub(jacobian);
    double errorSum = trueJacobian.sum();
    // std::cout << errorSum << "\n";
    BOOST_CHECK_SMALL(errorSum, epsilon);
  }
}

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorSurplusAdaptive) {
  std::vector<double> epsilons{1e-12};
  size_t dim = 3;
  size_t m = 2;
  size_t maxNumGridPoints = 100;
  size_t initialLevel = 0;
  size_t refinementsNum = 10;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  std::vector<size_t> degrees{3};
  for (size_t t = 0; t < degrees.size(); t++) {
    size_t degree = degrees[t];
    DataVector lb = testFunction->getLowerBounds();
    DataVector ub = testFunction->getUpperBounds();
    sgpp::optimization::SparseGridResponseSurfaceBsplineVector reSurf(testFunction, lb, ub,
                                                                      gridType, degree);
    reSurf.surplusAdaptive(maxNumGridPoints, initialLevel, refinementsNum);

    DataVector componentwiseErrors(m);
    size_t numMCPoints = 1000;
    DataVector errorVec = reSurf.averageL2Error(testFunction, componentwiseErrors, numMCPoints);
    double averageL2Err = errorVec[0];
    // std::cout << "l2 error: " << l2Error << "\n";
    // std::cout << "componentwise errors:\n" << componentwiseErrors.toString() << "\n";
    BOOST_CHECK_SMALL(averageL2Err, epsilons[t]);
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
    DataVector l2Error = reSurf.averageL2Error(testFunction, componentwiseErrors, numMCPoints);
    double averageL2Err = l2Error[0];
    // std::cout << "l2 error: " << l2Error << "\n";
    // std::cout << "componentwise errors:\n" << componentwiseErrors.toString() << "\n";
    BOOST_CHECK_SMALL(averageL2Err, epsilon);
  }
}

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorNRMSE) {
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

    DataMatrix componentwiseErrors(4, m);
    size_t numMCPoints = 1000;
    DataVector errorVec = reSurf.averageNRMSE(testFunction, componentwiseErrors, numMCPoints);
    double averageNRMSE = errorVec[0];
    double averageL2Err = errorVec[1];
    // std::cout << "NRMSE: " << averageNRMSE << "\n";
    // std::cout << "l2 error: " << averageL2Err << "\n";
    // std::cout << "componentwise errors:\n" << componentwiseErrors.toString() << "\n";
    BOOST_CHECK_SMALL(averageL2Err, epsilon);
    BOOST_CHECK_SMALL(averageNRMSE, epsilon);
  }
}

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorIntegral) {
  std::vector<double> epsilons{1e-12};
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  std::vector<size_t> degrees{3};
  for (size_t t = 0; t < degrees.size(); t++) {
    size_t degree = degrees[t];
    DataVector lb = testFunction->getLowerBounds();
    DataVector ub = testFunction->getUpperBounds();
    sgpp::optimization::SparseGridResponseSurfaceBsplineVector reSurf(testFunction, lb, ub,
                                                                      gridType, degree);
    reSurf.regular(level);
    sgpp::base::DataVector integrals = reSurf.getIntegrals();
    // std::cout << integrals.toString() << "\n";
    DataVector realIntegrals(2);
    realIntegrals[0] = 0.0;
    realIntegrals[1] = 3584.0 / 3.0;
    // std::cout << realIntegrals.toString() << "\n";
    realIntegrals.sub(integrals);
    // std::cout << realIntegrals.toString() << "\n";
    BOOST_CHECK_SMALL(fabs(realIntegrals[0]), epsilons[t]);
    BOOST_CHECK_SMALL(fabs(realIntegrals[1]), epsilons[t]);
  }
}

BOOST_AUTO_TEST_CASE(testResponseSurfaceBsplineVectorSerialize) {
  double epsilon = 1e-14;
  size_t dim = 3;
  size_t m = 2;
  size_t level = 1;  // small level to keep I/O operatiosn quick
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  size_t degree = 3;
  DataVector lb = testFunction->getLowerBounds();
  DataVector ub = testFunction->getUpperBounds();
  sgpp::optimization::SparseGridResponseSurfaceBsplineVector reSurf(testFunction, lb, ub, gridType,
                                                                    degree);
  reSurf.regular(level);

  // serialize
  std::string gridStr = reSurf.serializeGrid();
  DataMatrix coeffs = reSurf.getCoefficients();
  coeffs.toFile("testCoefffs.dat");
  std::ofstream out("testGrid.dat");
  out << gridStr;
  out.close();

  sgpp::optimization::SparseGridResponseSurfaceBsplineVector loadedReSurf(
      dim, m, lb, ub, "testGrid.dat", degree, "testCoefffs.dat");

  DataVector point(dim, 0.337);
  DataVector reSurfEval = reSurf.eval(point);
  DataVector loadedReSurfEval = loadedReSurf.eval(point);
  DataVector trueEval(m);
  testFunction->eval(point, trueEval);

  reSurfEval.sub(trueEval);
  loadedReSurfEval.sub(trueEval);

  double evalErr = reSurfEval.sum();
  double loadedEvalErr = loadedReSurfEval.sum();
  // std::cout << "evalErr = " << evalErr << "\n";
  BOOST_CHECK_SMALL(fabs(evalErr), epsilon);
  BOOST_CHECK_SMALL(fabs(loadedEvalErr), epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/DistributionLogNormal.hpp>
#include <sgpp/base/tools/DistributionNormal.hpp>
#include <sgpp/base/tools/DistributionTruncNormal.hpp>
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/function/scalar/SplineResponseSurface.hpp>
#include <sgpp/optimization/function/vector/SplineResponseSurfaceVector.hpp>

// debug
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::ScalarFunction;
using sgpp::base::VectorFunction;

// res = x_0^0 + x_1^1 + x_2^2 +...
class scalarTestFunction : public ScalarFunction {
 public:
  explicit scalarTestFunction(size_t d) : ScalarFunction(d) {}

  double eval(const DataVector& x) override {
    double res = 0;
    for (size_t i = 0; i < this->d; i++) {
      res += std::pow(x[i], i);
    }
    return res;
  }

  sgpp::base::DataVector evalGradient(sgpp::base::DataVector v) {
    sgpp::base::DataVector gradient(this->d);
    gradient[0] = 0;
    for (size_t i = 1; i < this->d; i++) {
      gradient[i] = static_cast<double>(i) * std::pow(v[i], i - 1);
    }
    return gradient;
  }

  virtual void clone(std::unique_ptr<sgpp::base::ScalarFunction>& clone) const override {
    clone = std::unique_ptr<sgpp::base::ScalarFunction>(new scalarTestFunction(*this));
  }
  DataVector getLowerBounds() { return DataVector(d, -2.0); }
  DataVector getUpperBounds() { return DataVector(d, 2.0); }
};

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

// Test the relevant supported spline basis types
std::vector<sgpp::base::GridType> getGridTypes() {
  std::vector<sgpp::base::GridType> gridTypes;
  gridTypes.push_back(sgpp::base::GridType::NakBsplineBoundary);
  // gridTypes.push_back(sgpp::base::GridType::ModNakBspline);
  gridTypes.push_back(sgpp::base::GridType::NakBsplineExtended);
  gridTypes.push_back(sgpp::base::GridType::NakPBspline);
  return gridTypes;
}

#ifdef USE_EIGEN
BOOST_AUTO_TEST_SUITE(TestResponseSurfaces)

/**
 * ***********************************
 * SplineResponseSurface tests
 * ***********************************
 */

// create regular SplineResponseSurface and evaluate
BOOST_AUTO_TEST_CASE(testRegularSplineResponseSurfaceEval) {
  // default values
  std::vector<double> epsilons{0.06, 1e-14, 1e-14};
  size_t dim = 3;
  size_t level = 3;
  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    // replace default values for specific cases
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilons.insert(epsilons.begin(), {0.06, 7e-4, 5e-5});
    }

    std::vector<size_t> degrees{1, 3, 5};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
      reSurf.regular(level);

      DataVector point(dim, 0.337);
      double reSurfEval = reSurf.eval(point);
      double trueEval = testFunction->eval(point);
      double diff = reSurfEval - trueEval;
      // std::cout << "evalErr = " << diff << "\n";
      BOOST_CHECK_SMALL(fabs(diff), epsilons[t]);
    }
  }
}

// create adaptive SplineResponseSurface and evaluate
BOOST_AUTO_TEST_CASE(testSurplusAdaptiveSplineResponseSurfaceEval) {
  std::vector<double> epsilons{0.3, 1e-14, 1e-14};
  size_t dim = 3;
  size_t maxNumGridPoints = 30;
  size_t initialLevel = 1;
  size_t refinementsNum = 1;
  bool verbose = false;
  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    // replace default values for specific cases
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilons.insert(epsilons.begin(), {0.3, 0.0007, 0.023});
    }
    std::vector<size_t> degrees{1, 3, 5};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
      reSurf.surplusAdaptive(maxNumGridPoints, initialLevel, refinementsNum, verbose);

      DataVector point(dim, 0.337);
      double reSurfEval = reSurf.eval(point);
      double trueEval = testFunction->eval(point);
      double diff = reSurfEval - trueEval;
      // std::cout << "evalErr = " << diff << "\n";
      BOOST_CHECK_SMALL(fabs(diff), epsilons[t]);
    }
  }
}

// create regular SplineResponseSurface and evaluate jacobian
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceEvalGradient) {
  double epsilon = 1e-13;
  size_t dim = 3;
  size_t level = 3;
  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilon = 0.003;
    }
    std::vector<size_t> degrees{3};
    for (auto& degree : degrees) {
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
      reSurf.regular(level);

      DataVector point(dim, 1.0 / 6.0);
      DataVector jacobian(dim);
      double reSurfEval = reSurf.evalGradient(point, jacobian);
      // suppress unused variable warning of gcc by doing something with reSurfEval
      reSurfEval += 1;
      DataVector trueJacobian = testFunction->evalGradient(point);
      // std::cout << trueJacobian.toString() << "\n";
      // std::cout << jacobian.toString() << "\n";
      trueJacobian.sub(jacobian);
      double errorSum = trueJacobian.sum();
      // std::cout << errorSum << "\n";
      BOOST_CHECK_SMALL(errorSum, epsilon);
    }
  }
}
// create regular SplineResponseSurface and integrate
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceIntegral) {
  double epsilon = 1e-12;
  size_t dim = 3;
  size_t level = 3;
  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilon = 1.2;
    }
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
      reSurf.regular(level);
      double integral = reSurf.getIntegral();
      // std::cout << integral << "\n";
      double realIntegral = 448.0 / 3.0;
      // std::cout << realIntegral << "\n";
      double diff = fabs(realIntegral - integral);
      // std::cout << diff << "\n";
      BOOST_CHECK_SMALL(diff, epsilon);
    }
  }
}

// create Ritter Novak SplineResponseSurface and find minimum
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceOptimization) {
  double epsilon = 1e-12;
  size_t dim = 3;
  // Ritter Novak parameters
  size_t maxNumGridPoints = 20;
  double gamma = 0.85;
  size_t initialLevel = 1;

  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();

  for (auto& gridType : gridTypes) {
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilon = 0.002;
    }
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
      reSurf.ritterNovak(maxNumGridPoints, gamma, initialLevel);
      DataVector argmin = reSurf.optimize();
      // This function has no unique minimum. [x,-2,0] is minimal for all x.
      // So we compare if the objective function evaluated at the found minimum is
      // close to the actual value
      // std::cout << argmin << "\n";
      double realMinimumValue = -1.0;
      double foundMinimumValue = testFunction->eval(argmin);
      // std::cout << foundMinimumValue << "\n";
      double diff = fabs(foundMinimumValue - realMinimumValue);
      // std::cout << diff << "\n";
      BOOST_CHECK_SMALL(diff, epsilon);
    }
  }
}

// create regular SplineResponseSurface and calculate mean
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceMean) {
  double epsilon = 1e-5;
  size_t dim = 3;
  size_t level = 3;
  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilon = 0.007;
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      epsilon = 0.002;
    } else if (gridType == sgpp::base::GridType::NakPBspline) {
      epsilon = 0.0017;
    }
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
      reSurf.regular(level);
      sgpp::base::DistributionsVector pdfs;
      auto normalpdf = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
      auto truncnormalpdf = std::make_shared<sgpp::base::DistributionTruncNormal>(0, 1, -2, 2);
      auto uniformpdf = std::make_shared<sgpp::base::DistributionUniform>(-2, 2);
      pdfs.push_back(uniformpdf);
      pdfs.push_back(truncnormalpdf);
      pdfs.push_back(truncnormalpdf);
      size_t quadOrder = 5;
      double mean = reSurf.getMean(pdfs, quadOrder);
      // std::cout << "mean: " << mean << " \n";
      // The reference mean is only determined numerically with nakBsplineBoundary level 7 (3809
      // grid points), trusting that the routine worked at this point in time For a foolproof test,
      // one might want to claulcate the mean by hand.
      double realMean = 1.61769;
      double diff = fabs(realMean - mean);
      // std::cout << diff << "\n";
      BOOST_CHECK_SMALL(diff, epsilon);
    }
  }
}

// ToDo(rehmemk) This must be debugged! Currently the variance is negativ!
// create regular SplineResponseSurface and calculate variance
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVariance) {
  double epsilon = 1e-5;
  size_t dim = 3;
  size_t level = 3;
  auto testFunction = std::make_shared<scalarTestFunction>(dim);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    if (gridType == sgpp::base::GridType::ModNakBspline) {
      epsilon = 0.05;
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      epsilon = 0.006;
    }
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurface reSurf(testFunction, lb, ub, gridType, degree);

      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";

      reSurf.regular(level);
      sgpp::base::DistributionsVector pdfs;
      auto truncnormalpdf = std::make_shared<sgpp::base::DistributionTruncNormal>(0, 1, -2, 2);
      auto uniformpdf = std::make_shared<sgpp::base::DistributionUniform>(-2, 2);
      pdfs.push_back(uniformpdf);
      pdfs.push_back(truncnormalpdf);
      pdfs.push_back(truncnormalpdf);
      size_t quadOrder = 5;
      DataVector varVec = reSurf.getVariance(pdfs, quadOrder);
      double var = varVec[0];
      double meanSquare = varVec[1];
      double mean = varVec[2];

      // std::cout << var << "\n";
      // std::cout << meanSquare << "\n";
      // std::cout << mean << "\n";

      // The reference mean is only determined numerically with nakBsplineBoundary level 7 (3809
      // grid points), trusting that the routine worked at this point in time For a foolproof test,
      // one might want to calculate the mean and variance by hand.
      double realVar = 1.70451;
      double realMeanSquare = 4.32144;
      double realMean = 1.61769;
      double diffVar = fabs(realVar - var);
      double diffMeanSquare = fabs(realMeanSquare - meanSquare);
      double diffMean = fabs(realMean - mean);

      // std::cout << diffMean << "\n";

      BOOST_CHECK_SMALL(diffVar, epsilon);
      BOOST_CHECK_SMALL(diffMeanSquare, epsilon);
      BOOST_CHECK_SMALL(diffMean, epsilon);
    }
  }
}

/**
 * ***********************************
 * SplineResponseSurfaceVector tests
 * ***********************************
 */

// create regular SplineResponseSurfaceVector and evaluate
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorEval) {
  std::vector<double> epsilons{0.8, 1e-13, 1e-13};
  size_t dim = 3;
  size_t m = 2;
  size_t level = 4;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{1, 3, 5};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
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
}

// create regular SplineResponseSurfaceVector and evaluate jacobian
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorEvalJacobian) {
  double epsilon = 1e-13;
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (auto& degree : degrees) {
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
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
}

// create adaptive SplineResponseSurfaceVector and call L2 error routine
BOOST_AUTO_TEST_CASE(testSurplusAdaptiveSplineResponseSurfaceVector) {
  std::vector<double> epsilons{1e-12};
  size_t dim = 3;
  size_t m = 2;
  size_t maxNumGridPoints = 100;
  size_t initialLevel = 1;
  size_t refinementsNum = 10;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
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
}

// create regular SplineResponseSurfaceVector and call L2 error routine
BOOST_AUTO_TEST_CASE(testRegularSplineResponseSurfaceVectorL2) {
  double epsilon = 1e-11;
  size_t dim = 3;
  size_t m = 3;
  size_t level = 4;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (auto& degree : degrees) {
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
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
}

// create regular SplineResponseSurfaceVector and call NRMSE error routine
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorNRMSE) {
  double epsilon = 1e-11;
  size_t dim = 3;
  size_t m = 3;
  size_t level = 4;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (auto& degree : degrees) {
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
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
}

// create regular SplineResponseSurfaceVector and integrate
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorIntegral) {
  std::vector<double> epsilons{1e-11};
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing " << reSurf.getGrid()->getTypeAsString() << " of deg" << degree <<
      // "\n";
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
}

// create regular SplineResponseSurfaceVector and calculate mean
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorMean) {
  std::vector<double> epsilons{1e-12};
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing testSplineResponseSurfaceVectorMean "
      //           << reSurf.getGrid()->getTypeAsString() << " of deg" << degree << "\n";
      reSurf.regular(level);

      // sgpp::base::DataMatrix componentwiseErrors;
      // size_t numMCPoints = 1000;
      // DataVector errorVec = reSurf.averageNRMSE(testFunction, componentwiseErrors, numMCPoints);
      // double averageNRMSE = errorVec[0];
      // double averageL2Err = errorVec[1];
      // std::cout << "average NRMSE: " << averageNRMSE << "   average L2:" << averageL2Err << "\n";

      sgpp::base::DistributionsVector pdfs(0);
      auto truncnormalpdf = std::make_shared<sgpp::base::DistributionTruncNormal>(0, 1, -2, 2);
      auto uniformpdf = std::make_shared<sgpp::base::DistributionUniform>(-2, 2);
      pdfs.push_back(truncnormalpdf);
      pdfs.push_back(truncnormalpdf);
      pdfs.push_back(uniformpdf);

      size_t quadOrder = 15;
      sgpp::base::DataVector means = reSurf.getMeans(pdfs, quadOrder);
      // std::cout << "means: " << means.toString() << "\n";
      // These are only determined numerically, trusting that the routine worked at this point in
      // time. (nakBsplineBoundary, level 7, quadrature order 100 -> NRMSE 6e-17).
      //  For a foolproof test, one might want to claulcate some means by hand.
      DataVector realMeans(2);
      realMeans[0] = 0.0;
      realMeans[1] = 1.44574984199963587628e+01;
      // std::cout << realMeans.toString() << "\n";
      realMeans.sub(means);
      // std::cout << "mean diffs: " << realMeans.toString() << "\n";
      BOOST_CHECK_SMALL(fabs(realMeans[0]), epsilons[t]);
      BOOST_CHECK_SMALL(fabs(realMeans[1]), epsilons[t]);
    }
  }
}

// create regular SplineResponseSurfaceVector and calculate variance
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorVariance) {
  std::vector<double> epsilons{1e-11};
  size_t dim = 3;
  size_t m = 2;
  size_t level = 3;
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  std::vector<sgpp::base::GridType> gridTypes = getGridTypes();
  for (auto& gridType : gridTypes) {
    std::vector<size_t> degrees{3};
    for (size_t t = 0; t < degrees.size(); t++) {
      size_t degree = degrees[t];
      DataVector lb = testFunction->getLowerBounds();
      DataVector ub = testFunction->getUpperBounds();
      sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType,
                                                             degree);
      // std::cout << "Testing testSplineResponseSurfaceVectorVariance "
      //           << reSurf.getGrid()->getTypeAsString() << " of deg" << degree << "\n";
      reSurf.regular(level);

      // sgpp::base::DataMatrix componentwiseErrors;
      // size_t numMCPoints = 1000;
      // DataVector errorVec = reSurf.averageNRMSE(testFunction, componentwiseErrors, numMCPoints);
      // double averageNRMSE = errorVec[0];
      // double averageL2Err = errorVec[1];
      // std::cout << "average NRMSE: " << averageNRMSE << "   average L2:" << averageL2Err << "\n";

      sgpp::base::DistributionsVector pdfs(0);
      auto truncnormalpdf = std::make_shared<sgpp::base::DistributionTruncNormal>(0, 1, -2, 2);
      auto uniformpdf = std::make_shared<sgpp::base::DistributionUniform>(-2, 2);
      pdfs.push_back(truncnormalpdf);
      pdfs.push_back(truncnormalpdf);
      pdfs.push_back(uniformpdf);
      size_t quadOrder = 15;
      sgpp::base::DataVector means(m);
      sgpp::base::DataVector meanSquares(m);
      sgpp::base::DataVector variances = reSurf.getVariances(pdfs, quadOrder, means, meanSquares);
      // std::cout << "variances: " << variances.toString() << "\n";

      // These are only determined numerically, trusting that the routine worked at this point in
      // time.
      // Calculated with NakBsplineBoundary level 7, quadOrder 100
      // For a foolproof test, one might want to claulcate some variances by hand.
      DataVector realVariances(m);
      realVariances[0] = 1.44574984199962521814e+01;
      realVariances[1] = 3.15930288844802475978e+02;
      // realVariances[2] = 2.65851725912918009271e+04;
      // std::cout << realVariances.toString() << "\n";
      realVariances.sub(variances);
      // std::cout << "variance diff: " << realVariances.toString() << "\n";
      BOOST_CHECK_SMALL(fabs(realVariances[0]), epsilons[t]);
      BOOST_CHECK_SMALL(fabs(realVariances[1]), epsilons[t]);
    }
  }
}

// create regular SplineResponseSurfaceVector, serialize it to files and create
// new  SplineResponseSurfaceVector by unserializing from files
BOOST_AUTO_TEST_CASE(testSplineResponseSurfaceVectorSerialize) {
  double epsilon = 1e-13;
  size_t dim = 3;
  size_t m = 2;
  size_t level = 1;  // small level to keep I/O operations quick
  auto testFunction = std::make_shared<multivariateTestFunction>(dim, m);
  sgpp::base::GridType gridType = sgpp::base::GridType::NakBsplineBoundary;
  size_t degree = 3;
  DataVector lb = testFunction->getLowerBounds();
  DataVector ub = testFunction->getUpperBounds();
  sgpp::optimization::SplineResponseSurfaceVector reSurf(testFunction, lb, ub, gridType, degree);
  reSurf.regular(level);

  // serialize
  std::string gridStr = reSurf.serializeGrid();
  DataMatrix coeffs = reSurf.getCoefficients();
  std::string coeffFileName = "reSurf_testCoeffs.dat";
  std::string gridFileName = "reSurf_testGrid.dat";
  coeffs.toFile(coeffFileName);
  std::ofstream out(gridFileName);
  out << gridStr;
  out.close();

  sgpp::optimization::SplineResponseSurfaceVector loadedReSurf(dim, m, lb, ub, gridFileName, degree,
                                                               coeffFileName);

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
  remove("reSurf_testCoeffs.dat");
  remove("reSurf_testGrid.dat");
}

BOOST_AUTO_TEST_SUITE_END()
#endif

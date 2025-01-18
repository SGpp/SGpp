// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp>
#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaOptimization.hpp>
#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaTransformation.hpp>
#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaVertexMethod.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>
#include <sgpp/optimization/fuzzy/FuzzyIntervalViaConfidenceInterval.hpp>
#include <sgpp/optimization/fuzzy/FuzzyIntervalViaMembershipFunction.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>
#include <sgpp/optimization/fuzzy/QuasiGaussianFuzzyNumber.hpp>
#include <sgpp/optimization/fuzzy/TriangularFuzzyInterval.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>

#include <vector>

using sgpp::optimization::FuzzyExtensionPrinciple;
using sgpp::optimization::FuzzyExtensionPrincipleViaOptimization;
using sgpp::optimization::FuzzyExtensionPrincipleViaTransformation;
using sgpp::optimization::FuzzyExtensionPrincipleViaVertexMethod;
using sgpp::optimization::FuzzyInterval;
using sgpp::optimization::FuzzyIntervalViaConfidenceInterval;
using sgpp::optimization::FuzzyIntervalViaMembershipFunction;
using sgpp::optimization::InterpolatedFuzzyInterval;
using sgpp::optimization::QuasiGaussianFuzzyNumber;
using sgpp::optimization::TriangularFuzzyInterval;

class TestFuzzyInterval1 : public FuzzyIntervalViaMembershipFunction {
 public:
  TestFuzzyInterval1() : FuzzyIntervalViaMembershipFunction(-2.0, 4.0, 1.0, 3.0) {}

  double evaluateMembershipFunction(double x) const override {
    if (x < -2.0) {
      return 0.0;
    } else if (x < 1.0) {
      return (x + 2.0) / 3.0;
    } else if (x < 3.0) {
      return 1.0;
    } else if (x < 4.0) {
      return 4.0 - x;
    } else {
      return 0.0;
    }
  }
};

class TestFuzzyInterval2 : public FuzzyIntervalViaConfidenceInterval {
 public:
  TestFuzzyInterval2() : FuzzyIntervalViaConfidenceInterval(-2.0, 4.0) {}

  double evaluateConfidenceIntervalLowerBound(double alpha) const override {
    return 3.0 * alpha - 2.0;
  }

  double evaluateConfidenceIntervalUpperBound(double alpha) const override {
    return 4.0 - alpha;
  }
};

class BilinearFunction : public sgpp::base::ScalarFunction {
 public:
  BilinearFunction() : ScalarFunction(2) {}
  ~BilinearFunction() override {}

  inline double eval(const sgpp::base::DataVector& x) override {
    return (8.0 * x[0]) * (8.0 * x[1]) / 10.0;
  }

  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(new BilinearFunction());
  }
};

class BilinearFunctionGradient : public sgpp::base::ScalarFunctionGradient {
 public:
  BilinearFunctionGradient() : ScalarFunctionGradient(2) {}
  ~BilinearFunctionGradient() override {}

  inline double eval(const sgpp::base::DataVector& x, sgpp::base::DataVector& gradient) override {
    gradient.resize(2);
    gradient[0] = 8.0 * (8.0 * x[1]) / 10.0;
    gradient[1] = 8.0 * (8.0 * x[0]) / 10.0;
    return (8.0 * x[0]) * (8.0 * x[1]) / 10.0;
  }

  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionGradient>(new BilinearFunctionGradient());
  }
};

BOOST_AUTO_TEST_CASE(TestFuzzyIntervalBinarySearch) {
  TestFuzzyInterval1 interval1;
  TestFuzzyInterval2 interval2;

  interval1.setBinarySearchTolerance(1e-10);
  BOOST_CHECK_EQUAL(interval1.getBinarySearchTolerance(), 1e-10);
  interval2.setBinarySearchTolerance(1e-10);
  BOOST_CHECK_EQUAL(interval2.getBinarySearchTolerance(), 1e-10);

  for (size_t i = 0; i <= 20; i++) {
    const double x = 7.0 * static_cast<double>(i) / 20.0 - 3.0;
    BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(x),
                      interval2.evaluateMembershipFunction(x), 1e-6);
  }

  for (size_t i = 0; i <= 20; i++) {
    const double alpha = static_cast<double>(i) / 20.0;
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(alpha),
                      interval2.evaluateConfidenceIntervalLowerBound(alpha), 1e-6);
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(alpha),
                      interval2.evaluateConfidenceIntervalUpperBound(alpha), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(TestTriangularFuzzyInterval) {
  {
    const TriangularFuzzyInterval interval(1.2, 3.4);
    BOOST_CHECK_CLOSE(interval.getSupportLowerBound(), -2.2, 1e-8);
    BOOST_CHECK_CLOSE(interval.getSupportUpperBound(), 4.6, 1e-8);
    BOOST_CHECK_EQUAL(interval.getRightMean(), 1.2);
    BOOST_CHECK_EQUAL(interval.getLeftMean(), 1.2);
    BOOST_CHECK_EQUAL(interval.getRightMean(), 1.2);
    BOOST_CHECK_EQUAL(interval.getLeftSpread(), 3.4);
    BOOST_CHECK_EQUAL(interval.getRightSpread(), 3.4);
  }

  {
    const TriangularFuzzyInterval interval(1.2, 3.4, 5.6);
    BOOST_CHECK_CLOSE(interval.getSupportLowerBound(), -2.2, 1e-8);
    BOOST_CHECK_CLOSE(interval.getSupportUpperBound(), 6.8, 1e-8);
    BOOST_CHECK_EQUAL(interval.getLeftMean(), 1.2);
    BOOST_CHECK_EQUAL(interval.getRightMean(), 1.2);
    BOOST_CHECK_EQUAL(interval.getLeftSpread(), 3.4);
    BOOST_CHECK_EQUAL(interval.getRightSpread(), 5.6);
  }

  {
    TriangularFuzzyInterval interval1(1.2, 3.4, 5.6, 7.8);
    TriangularFuzzyInterval interval2(interval1);

    for (const TriangularFuzzyInterval* interval : {&interval1, &interval2}) {
      BOOST_CHECK_CLOSE(interval->getSupportLowerBound(), -4.4, 1e-8);
      BOOST_CHECK_CLOSE(interval->getSupportUpperBound(), 11.2, 1e-8);
      BOOST_CHECK_EQUAL(interval->getLeftMean(), 1.2);
      BOOST_CHECK_EQUAL(interval->getRightMean(), 3.4);
      BOOST_CHECK_EQUAL(interval->getLeftSpread(), 5.6);
      BOOST_CHECK_EQUAL(interval->getRightSpread(), 7.8);
    }

    BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(1.2), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(3.4), 1.0, 1e-8);
    BOOST_CHECK_SMALL(interval1.evaluateMembershipFunction(-4.4), 1e-8);
    BOOST_CHECK_SMALL(interval1.evaluateMembershipFunction(11.2), 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(-3.7), 0.125, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(7.3), 0.5, 1e-8);

    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(0.0), -4.4, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(0.0), 11.2, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(1.0), 1.2, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(1.0), 3.4, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(0.5), -1.6, 1e-8);
    BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(0.5), 7.3, 1e-8);

    interval1.setNumberOfIntegralSamples(12345);
    BOOST_CHECK_EQUAL(interval1.getNumberOfIntegralSamples(), 12345);

    BOOST_CHECK_CLOSE(
        interval1.computeL1Norm(TriangularFuzzyInterval::NormMode::ViaMembershipFunction),
        8.9, 1e-2);
    BOOST_CHECK_CLOSE(
        interval1.computeL1Norm(TriangularFuzzyInterval::NormMode::ViaConfidenceInterval),
        8.9, 1e-2);
    BOOST_CHECK_CLOSE(
        interval1.computeL2Norm(TriangularFuzzyInterval::NormMode::ViaMembershipFunction),
        2.58198890, 1e-2);
    BOOST_CHECK_CLOSE(
        interval1.computeL2Norm(TriangularFuzzyInterval::NormMode::ViaConfidenceInterval),
        9.70429458, 1e-2);
    BOOST_CHECK_CLOSE(
        interval1.computeLinfNorm(TriangularFuzzyInterval::NormMode::ViaMembershipFunction),
        1.0, 1e-2);
    BOOST_CHECK_CLOSE(
        interval1.computeLinfNorm(TriangularFuzzyInterval::NormMode::ViaConfidenceInterval),
        15.6, 1e-2);

    for (TriangularFuzzyInterval::NormMode normMode : {
        TriangularFuzzyInterval::NormMode::ViaMembershipFunction,
        TriangularFuzzyInterval::NormMode::ViaConfidenceInterval,
    }) {
      BOOST_CHECK_SMALL(interval1.computeL1Error(interval2, normMode), 1e-6);
      BOOST_CHECK_SMALL(interval1.computeL2Error(interval2, normMode), 1e-6);
      BOOST_CHECK_SMALL(interval1.computeLinfError(interval2, normMode), 1e-6);
      BOOST_CHECK_SMALL(interval1.computeRelativeL1Error(interval2, normMode), 1e-6);
      BOOST_CHECK_SMALL(interval1.computeRelativeL2Error(interval2, normMode), 1e-6);
      BOOST_CHECK_SMALL(interval1.computeRelativeLinfError(interval2, normMode), 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestQuasiGaussianFuzzyNumber) {
  const QuasiGaussianFuzzyNumber number1(1.2, 3.4, 0.5);
  const QuasiGaussianFuzzyNumber number2(number1);

  for (const QuasiGaussianFuzzyNumber* number : {&number1, &number2}) {
    BOOST_CHECK_CLOSE(number->getSupportLowerBound(), -0.5, 1e-8);
    BOOST_CHECK_CLOSE(number->getSupportUpperBound(), 2.9, 1e-8);
    BOOST_CHECK_EQUAL(number->getMean(), 1.2);
    BOOST_CHECK_EQUAL(number->getStdev(), 3.4);
    BOOST_CHECK_EQUAL(number->getCutoff(), 0.5);
  }

  BOOST_CHECK_SMALL(number1.evaluateMembershipFunction(-0.50001), 1e-8);
  BOOST_CHECK_SMALL(number1.evaluateMembershipFunction(2.90001), 1e-8);
  BOOST_CHECK_CLOSE(number1.evaluateMembershipFunction(1.2), 1.0, 1e-8);
  BOOST_CHECK_CLOSE(number1.evaluateMembershipFunction(-0.49999), 0.882498, 1e-4);
  BOOST_CHECK_CLOSE(number1.evaluateMembershipFunction(2.89999), 0.882498, 1e-4);

  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalLowerBound(0.0), -0.5, 1e-8);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalUpperBound(0.0), 2.9, 1e-8);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalLowerBound(0.3), -0.5, 1e-4);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalUpperBound(0.3), 2.9, 1e-4);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalLowerBound(0.9), -0.360748, 1e-4);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalUpperBound(0.9), 2.76075, 1e-4);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalLowerBound(1.0), 1.2, 1e-8);
  BOOST_CHECK_CLOSE(number1.evaluateConfidenceIntervalUpperBound(1.0), 1.2, 1e-8);
}

BOOST_AUTO_TEST_CASE(TestInterpolatedFuzzyInterval) {
  const sgpp::base::DataVector xData(std::vector<double>({1.2, 3.4, 5.6, 7.8, 9.10, 11.12}));
  const sgpp::base::DataVector alphaData(std::vector<double>({0.0, 0.5, 1.0, 1.0, 0.25, 0.0}));

  BOOST_CHECK_CLOSE(InterpolatedFuzzyInterval::getCoreLowerBound(xData, alphaData), 5.6, 1e-4);
  BOOST_CHECK_CLOSE(InterpolatedFuzzyInterval::getCoreUpperBound(xData, alphaData), 7.8, 1e-4);

  InterpolatedFuzzyInterval interval1(xData, alphaData);
  InterpolatedFuzzyInterval interval2(interval1);
  TriangularFuzzyInterval interval3(1.2, 3.4);

  BOOST_CHECK(InterpolatedFuzzyInterval::tryDowncast(interval1) != nullptr);
  BOOST_CHECK(InterpolatedFuzzyInterval::tryDowncast(interval2) != nullptr);
  BOOST_CHECK(InterpolatedFuzzyInterval::tryDowncast(interval3) == nullptr);

  for (InterpolatedFuzzyInterval* interval : {&interval1, &interval2}) {
    const sgpp::base::DataVector& xData1(interval->getXData());
    const sgpp::base::DataVector& alphaData1(interval->getAlphaData());

    BOOST_CHECK_EQUAL(xData.size(), xData1.size());
    BOOST_CHECK_EQUAL(alphaData.size(), alphaData1.size());

    for (size_t i = 0; i < xData.size(); i++) {
      BOOST_CHECK_EQUAL(xData[i], xData1[i]);
      BOOST_CHECK_EQUAL(alphaData[i], alphaData1[i]);
    }
  }

  BOOST_CHECK_EQUAL(interval1.evaluateMembershipFunction(1.1), 0.0);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(3.4), 0.5, 1e-6);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(4.5), 0.75, 1e-6);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(6.0), 1.0, 1e-6);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(9.0), 0.30769231, 1e-6);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(10.11), 0.125, 1e-6);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(11.12), 0.0, 1e-6);
  BOOST_CHECK_CLOSE(interval1.evaluateMembershipFunction(12.0), 0.0, 1e-6);

  BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(0.0), 1.2, 1e-5);
  BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(0.0), 11.12, 1e-5);
  BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(1.0), 5.6, 1e-5);
  BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(1.0), 7.8, 1e-5);
  BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalLowerBound(0.5), 3.4, 1e-5);
  BOOST_CHECK_CLOSE(interval1.evaluateConfidenceIntervalUpperBound(0.5), 8.6666667, 1e-5);
}

BOOST_AUTO_TEST_CASE(TestFuzzyExtensionPrinciple) {
  sgpp::base::Printer::getInstance().setVerbosity(-1);
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(42);

  const BilinearFunction f;
  const BilinearFunctionGradient fGradient;
  sgpp::optimization::optimizer::AdaptiveGradientDescent optimizer(f, fGradient);

  std::vector<std::unique_ptr<FuzzyExtensionPrinciple>> principles;
  principles.push_back(std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaOptimization(f, 12)));
  principles.push_back(std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaOptimization(optimizer, 12)));
  principles.push_back(std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaTransformation(f, 12)));
  principles.push_back(std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaVertexMethod(f, 12)));

  TriangularFuzzyInterval xFuzzy1(0.625, 0.75, 0.25, 0.125);
  TriangularFuzzyInterval xFuzzy2(0.45, 0.55, 0.4, 0.4);
  std::vector<FuzzyInterval*> xFuzzy = {&xFuzzy1, &xFuzzy2};

  for (auto& principle : principles) {
    // test clone
    std::unique_ptr<FuzzyExtensionPrinciple> curPrinciple;
    principle->clone(curPrinciple);

    BOOST_CHECK_EQUAL(curPrinciple->getNumberOfAlphaSegments(), 12);
    const size_t m = 20;
    curPrinciple->setNumberOfAlphaSegments(m);
    BOOST_CHECK_EQUAL(curPrinciple->getNumberOfAlphaSegments(), m);

    std::unique_ptr<FuzzyInterval> yFuzzy(curPrinciple->apply(xFuzzy));

    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalLowerBound(0.0), 0.12, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalUpperBound(0.0),  5.32, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalLowerBound(0.25), 0.42, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalUpperBound(0.25), 4.59, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalLowerBound(0.5), 0.8, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalUpperBound(0.5), 3.9, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalLowerBound(0.75), 1.26, 5e0);
    BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalUpperBound(0.75), 3.25, 5e0);

    // the transformation method only places one point on the optimization domain for alpha = 1,
    // which will be inaccurate if the fuzzy input intervals are not fuzzy numbers, like here
    // ==> skip the transformation method for alpha = 1
    if (dynamic_cast<FuzzyExtensionPrincipleViaTransformation*>(principle.get()) == nullptr) {
      BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalLowerBound(1.0), 1.8, 1e0);
      BOOST_CHECK_CLOSE(yFuzzy->evaluateConfidenceIntervalUpperBound(1.0), 2.64, 1e0);
    }

    const sgpp::base::DataVector& alphaLevels = curPrinciple->getAlphaLevels();
    BOOST_CHECK_EQUAL(alphaLevels.size(), m + 1);
    BOOST_CHECK_EQUAL(alphaLevels[m/2], 0.5);

    const std::vector<sgpp::base::DataVector>& lowerBounds =
        curPrinciple->getOptimizationDomainsLowerBounds();
    BOOST_CHECK_EQUAL(lowerBounds.size(), m + 1);
    BOOST_CHECK_EQUAL(lowerBounds[m/2].size(), 2);
    BOOST_CHECK_CLOSE(lowerBounds[m/2][0], 0.5, 1e-4);
    BOOST_CHECK_CLOSE(lowerBounds[m/2][1], 0.25, 1e-4);

    const std::vector<sgpp::base::DataVector>& upperBounds =
        curPrinciple->getOptimizationDomainsUpperBounds();
    BOOST_CHECK_EQUAL(upperBounds.size(), m + 1);
    BOOST_CHECK_EQUAL(upperBounds[m/2].size(), 2);
    BOOST_CHECK_CLOSE(upperBounds[m/2][0], 0.8125, 1e-4);
    BOOST_CHECK_CLOSE(upperBounds[m/2][1], 0.75, 1e-4);

    const std::vector<sgpp::base::DataVector>& minimumPoints = curPrinciple->getMinimumPoints();
    BOOST_CHECK_EQUAL(minimumPoints.size(), m + 1);
    BOOST_CHECK_EQUAL(minimumPoints[m/2].size(), 2);
    BOOST_CHECK_CLOSE(minimumPoints[m/2][0], 0.5, 5e0);
    BOOST_CHECK_CLOSE(minimumPoints[m/2][1], 0.25, 5e0);

    const sgpp::base::DataVector& minimumValues = curPrinciple->getMinimumValues();
    BOOST_CHECK_EQUAL(minimumValues.size(), m + 1);
    BOOST_CHECK_CLOSE(minimumValues[m/2], 0.8, 5e0);

    const std::vector<sgpp::base::DataVector>& maximumPoints = curPrinciple->getMaximumPoints();
    BOOST_CHECK_EQUAL(maximumPoints.size(), m + 1);
    BOOST_CHECK_EQUAL(maximumPoints[m/2].size(), 2);
    BOOST_CHECK_CLOSE(maximumPoints[m/2][0], 0.8125, 5e0);
    BOOST_CHECK_CLOSE(maximumPoints[m/2][1], 0.75, 5e0);

    const sgpp::base::DataVector& maximumValues = curPrinciple->getMaximumValues();
    BOOST_CHECK_EQUAL(maximumValues.size(), m + 1);
    BOOST_CHECK_CLOSE(maximumValues[m/2], 3.9, 5e0);
  }
}

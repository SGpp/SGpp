// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGPPCombigridModule

#include <boost/test/unit_test.hpp>

#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

using sgpp::base::DataVector;
using sgpp::combigrid::AbstractMultiStorage;
using sgpp::combigrid::FloatArrayVector;
using sgpp::combigrid::CombigridMultiOperation;
using sgpp::combigrid::MultiFunction;
using sgpp::combigrid::Stopwatch;
using sgpp::combigrid::MCIntegrator;

double testFunction(DataVector const &coordinates) {
  double prod = 1.0;

  for (size_t i = 0; i < coordinates.getSize(); ++i) {
    double x = coordinates[i];
    prod *= std::exp(-x * x / static_cast<double>((i + 1) * (i + 1)));
  }

  return prod;
}

double testFunction2(DataVector const &coordinates) {
  double prod = 1.0;

  for (size_t i = 0; i < coordinates.getSize(); ++i) {
    double x = coordinates[i];
    prod *= std::cos(-x * x / static_cast<double>((i + 1) * (i + 1)));
  }

  return prod;
}

double testFunction3(DataVector const &coordinates) {
  double prod = 1.0;

  for (size_t i = 0; i < coordinates.getSize(); ++i) {
    double x = coordinates[i];
    prod *= (x * x) / (1 + x * x);
  }

  return prod;
}

double testFunction4(DataVector const &coordinates) {
  double prod = 1.0;

  for (size_t i = 0; i < coordinates.getSize(); ++i) {
    double x = coordinates[i];
    prod *= exp(-x * x);
  }

  return prod;
}

double testFunction5(DataVector const &x) {
  std::vector<double> k;
  for (size_t i = 0; i < x.getSize(); ++i) {
    double root = std::sqrt(static_cast<double>(i) +
                            6405.0);  // 6400 ist Wurzel von 80, die naechste Quadratzahl ist
                                      // 6561, also 161 mal Nachkommastellen ungleich 0
    double decimals = root - std::floor(root);
    k.push_back(static_cast<double>(i % 2 == 0 ? 1 : -1) * 100 * decimals *
                std::tan(static_cast<double>(i)));
  }

  double sum = 0.0;

  for (size_t i = 0; i < x.getSize(); ++i) {
    double tmp = 1;
    for (size_t j = 0; j < k.size(); ++j) {
      tmp *= std::sin((-1) * k[j] * (1 - x[i]) * static_cast<double>(i)) -
             sgpp::combigrid::pow(cos(k[k.size() - 1 - i] * x[i] * k[j]), 3) +
             std::sin((1 - x[i]) * (1 - x[i]) + x[x.getSize() - 1 - i]);
    }
    sum += tmp;
  }

  return sum;
}

double testFunction6(DataVector const &x) { return 1; }

double testFunction7(DataVector const &x) { return x[0]; }

double testFunctionAtan(DataVector const &x) {
  return std::atan(50 * (x[0] - .35)) + M_PI / 2 +
         4 * std::pow(x[1], 2);  // + exp(x[0] * x[1] - 1);
}

double testFunction8(DataVector const &coordinates) {
  double ans = 0.0;

  for (size_t i = 0; i < coordinates.getSize(); ++i) {
    ans += coordinates[i];
  }

  return ans;
}

/*
 void printCTResults(size_t d, size_t q) {
 const size_t samples = 10;
 auto ctInterpolator = CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(d,
 testFunction);
 auto domain = std::vector<std::pair<double, double>>(d, std::pair<double, double>(0.0,
 1.0));

 MCIntegrator integrator([&](DataVector const &x) -> double {
 double diff = testFunction(x) - ctInterpolator->evaluate(q, x);
 return diff * diff;
 });

 std::cout << "d = " << d << ", q = " << q << ": " << std::sqrt(integrator.average(domain,
 samples)) << std::endl;
 }
 */

void printDifferences(size_t d, std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> storage) {
  auto it = storage->getStoredDataIterator();

  std::cout << "Differences: \n";
  while (it->isValid()) {
    std::cout << "Level (";
    for (size_t i = 0; i < d - 1; ++i) {
      std::cout << it->indexAt(i);
      std::cout << ", ";
    }
    std::cout << it->indexAt(d - 1) << "): ";

    std::cout << it->value().norm() << "\n";

    it->moveToNext();
  }
}

void computeL2Error(size_t d, size_t q, std::shared_ptr<CombigridMultiOperation> ctInterpolator,
                    MultiFunction &func) {
  const size_t samples = 100;
  auto domain = std::vector<std::pair<double, double>>(d, std::pair<double, double>(0.0, 1.0));

  MCIntegrator integrator(
      [&ctInterpolator, q, func](std::vector<DataVector> const &params) -> DataVector {
        auto result = ctInterpolator->evaluate(q, params);
        // auto result = ctInterpolator->evaluateAdaptive(q * number, params);

        // printDifferences(d, ctInterpolator->getDifferences());

        for (size_t i = 0; i < params.size(); ++i) {
          double diff = func(params[i]) - result[i];
          result[i] = diff * diff;
        }

        return result;
      });

  std::cout << "d = " << d << ", q = " << q
            << ": err_l2 = " << std::sqrt(integrator.average(domain, samples)) << std::endl;
}

void computeQuadratureError(size_t d, size_t q,
                            std::shared_ptr<CombigridMultiOperation> ctInterpolator,
                            MultiFunction &func) {
  const size_t numSamples = 1e4;
  auto domain = std::vector<std::pair<double, double>>(d, std::pair<double, double>(0.0, 1.0));
  static std::default_random_engine generator(
      std::mt19937_64::default_seed);  // TODO(holzmudd): deterministic

  double result = 0.0;
  sgpp::base::DataVector coordinates(domain.size());
  for (size_t i = 0; i < numSamples; ++i) {
    for (size_t dim = 0; dim < domain.size(); ++dim) {
      std::uniform_real_distribution<double> distribution(domain[dim].first, domain[dim].second);
      coordinates[dim] = distribution(generator);
    }
    result += func(coordinates);
  }

  result /= static_cast<double>(numSamples);

  std::cout << "d = " << d << ", q = " << q
            << ": err_quad = " << std::abs(result - ctInterpolator->evaluate(q)[0]) << std::endl;
}

BOOST_AUTO_TEST_SUITE(testInterpolation)

BOOST_AUTO_TEST_CASE(testLinearInterpolation) {
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "Linear Interpolation" << std::endl;
  auto func = MultiFunction(testFunction3);
  for (size_t d = 2; d <= 5; ++d) {
    auto ctInterpolator =
        CombigridMultiOperation::createExpUniformBoundaryLinearInterpolation(d, func);
    std::cout << "- - - - - - - - - - - - - - " << std::endl;
    for (size_t w = 2; w <= 8; ++w) {
      Stopwatch stopwatch;
      computeL2Error(d, w, ctInterpolator, func);
      // stopwatch.log();
    }
  }
  auto quadrature = CombigridMultiOperation::createLinearLejaQuadrature(
      3, MultiFunction([](sgpp::base::DataVector const &x) { return 1.0; }));
  std::vector<DataVector> input(1, DataVector(0));
  auto result = quadrature->evaluate(3, input);
  std::cout << "Quadrature result: " << result[0] << "\n";
}

BOOST_AUTO_TEST_CASE(testExpChebyshevPolynomialInterpolation) {
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "Chebyshev Polynomial Interpolation" << std::endl;
  auto func = MultiFunction(testFunction3);
  for (size_t d = 2; d <= 5; ++d) {
    auto ctInterpolator =
        CombigridMultiOperation::createExpChebyshevPolynomialInterpolation(d, func);
    std::cout << "- - - - - - - - - - - - - - " << std::endl;
    for (size_t w = 2; w <= 8; ++w) {
      Stopwatch stopwatch;
      computeL2Error(d, w, ctInterpolator, func);
      // stopwatch.log();
    }
  }

  auto quadrature = CombigridMultiOperation::createLinearLejaQuadrature(
      3, MultiFunction([](sgpp::base::DataVector const &x) { return 1.0; }));
  std::vector<DataVector> input(1, DataVector(0));
  auto result = quadrature->evaluate(3, input);
  std::cout << "Quadrature result: " << result[0] << "\n";
}

BOOST_AUTO_TEST_CASE(testPolynomialInterpolation) {
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "Polynomial Interpolation" << std::endl;
  auto func = MultiFunction(testFunctionAtan);
  for (size_t d = 2; d <= 5; ++d) {
    std::cout << "- - - - - - - - - - - - - - " << std::endl;
    auto ctInterpolator =
        CombigridMultiOperation::createExpClenshawCurtisPolynomialInterpolation(d, func);
    for (size_t w = 2; w <= 8; ++w) {
      computeL2Error(d, w, ctInterpolator, func);
    }
    auto ctQuadrature = CombigridMultiOperation::createExpClenshawCurtisQuadrature(d, func);
    for (size_t w = 2; w <= 8; ++w) {
      computeQuadratureError(d, w, ctQuadrature, func);
    }
  }

  auto quadrature = CombigridMultiOperation::createLinearLejaQuadrature(
      3, MultiFunction([](sgpp::base::DataVector const &x) { return 1.0; }));
  std::vector<DataVector> input(1, DataVector(0));
  auto result = quadrature->evaluate(3, input);
  std::cout << "Quadrature result: " << result[0] << "\n";
}

BOOST_AUTO_TEST_CASE(testLinearL2LejaPolynomialInterpolation) {
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "L2 Leja Polynomial Interpolation" << std::endl;
  auto func = MultiFunction(testFunction3);
  for (size_t d = 2; d <= 5; ++d) {
    std::cout << "- - - - - - - - - - - - - - " << std::endl;
    auto ctInterpolator =
        CombigridMultiOperation::createLinearL2LejaPolynomialInterpolation(d, func, 2);
    for (size_t w = 2; w <= 8; ++w) {
      computeL2Error(d, w, ctInterpolator, func);
    }
    auto ctQuadrature = CombigridMultiOperation::createLinearL2LejaQuadrature(d, func, 2);
    for (size_t w = 2; w <= 8; ++w) {
      computeQuadratureError(d, w, ctQuadrature, func);
    }
  }
}

BOOST_AUTO_TEST_CASE(testBsplinedeg3Interpolation) {
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << "B-Spline Interpolation, degree=3\n" << std::endl;
  auto func = MultiFunction(testFunction3);
  for (size_t d = 2; d <= 5; ++d) {
    std::cout << "- - - - - - - - - - - - - - " << std::endl;
    size_t degree = 3;
    auto ctInterpolator =
        CombigridMultiOperation::createExpUniformBoundaryBsplineInterpolation(d, func, degree);
    for (size_t w = 2; w <= 6; ++w) {
      computeL2Error(d, w, ctInterpolator, func);
    }

    auto ctQuadrature =
        CombigridMultiOperation::createExpUniformBoundaryBsplineQuadrature(d, func, degree);
    for (size_t w = 2; w <= 6; ++w) {
      computeQuadratureError(d, w, ctQuadrature, func);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

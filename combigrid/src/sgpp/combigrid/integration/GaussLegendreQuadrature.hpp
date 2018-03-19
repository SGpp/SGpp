// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <algorithm>
#include <cmath>
#include <vector>
// #include <future>
#include <iostream>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace combigrid {

/**
 * Uses Gauss-Legendre quadrature to approximate the integral of a function on the domain [0, 1].
 */
class GaussLegendreQuadrature {
  base::DataVector roots;
  base::DataVector weights;

 public:
  /**
   * Constructor.
   * @param numPoints The number of grid points that should be used. This number should be > 0. If
   * not enough grid points are available, the maximum number of available grid points is taken.
   */
  explicit GaussLegendreQuadrature(size_t numPoints = 0);

  void initialize(size_t numPoints);

  template <typename Func>
  double evaluate(Func const &func, double a = 0.0, double b = 1.0) {
    double width = b - a;
    double sum = 0.0;
    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x_unit = roots[i], x_prob = a + width * roots[i];
      sum += weights[i] * func(x_unit, x_prob);
    }
    return width * sum;
  }

  //  template <typename Func>
  //  static double evaluate_parallel(Func const &func, double a = 0.0, double b = 1.0,
  //                                  size_t numGaussPoints = 1, size_t incrementQuadraturePoints =
  //                                  1,
  //                                  double tol = 1e-14, size_t numThreads = 0) {
  //    if (incrementQuadraturePoints > 1) {
  //      if (numThreads == 0) {
  // #ifdef _OPENMP
  //        numThreads = std::max(1, omp_get_max_threads());
  // #else
  //        numThreads = 1;
  // #endif
  //      }
  //      if (numThreads > 1) {
  //        double segmentWidth = (b - a) / static_cast<double>(numThreads);
  //
  //        std::vector<std::future<double>> futures;
  //        for (double xlower = a; xlower < b; xlower += segmentWidth) {
  //          futures.push_back(std::async(std::launch::async, evaluate_iteratively<Func>, func,
  //          xlower,
  //                                       xlower + segmentWidth, numGaussPoints,
  //                                       incrementQuadraturePoints, tol));
  //        }
  //        double sum = 0.0;
  //        for (auto &handle : futures) {
  //          sum += handle.get();
  //        }
  //        return sum;
  //      }
  //    }
  //
  //    return evaluate_iteratively(func, a, b, numGaussPoints, incrementQuadraturePoints, tol);
  //  }

  template <typename Func>
  double evaluate_iteratively(Func const &func, double a = 0.0, double b = 1.0,
                              size_t numGaussPoints = 0, size_t incrementQuadraturePoints = 1,
                              double tol = 1e-14) {
    double sum = 0.0, sum_old = 0.0;
    double width = b - a;

    auto &quadRule = base::GaussLegendreQuadRule1D::getInstance();

    // performing Gauss-Legendre integration
    double err = 1e14;
    size_t iteration = 0;

    // make sure that the quadrature is done at least once
    numGaussPoints = std::min(quadRule.getMaxSupportedLevel(), numGaussPoints);

    while (err > tol && numGaussPoints <= quadRule.getMaxSupportedLevel()) {
      quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, weights);

      sum = 0.0;
      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x_unit = roots[i], w = weights[i];
        double x_trans = a + width * x_unit;
        sum += w * func(x_unit, x_trans);
      }
      sum *= width;

      if (iteration > 0) {
        err = std::fabs(sum_old - sum);
      }

      sum_old = sum;
      numGaussPoints += incrementQuadraturePoints;
      iteration += 1;

      if (incrementQuadraturePoints == 0) {
        return sum;
      }
    }
    return sum;
  }

  template <typename Func>
  double evaluate_logsumexp(Func const &logfunc, SingleFunction weightFunction, double a = 0.0,
                            double b = 1.0) {
    return std::log(evaluate_sumexp(logfunc, weightFunction, a, b));
  }

  template <typename Func>
  double evaluate_sumexp(Func const &logfunc, SingleFunction weightFunction, double a = 0.0,
                         double b = 1.0) {
    double width = b - a;
    double logIntegral = 0.0;

    // performing Gauss-Legendre integration
    sgpp::base::DataVector psi(roots.size());
    // --------------------------------------------------------------
    // do quadrature
    double max_psi = -std::numeric_limits<double>::max();

    for (size_t i = 0; i < roots.getSize(); ++i) {
      double x_trans = a + width * roots[i];
      psi[i] = std::log(weights[i]) + std::log(weightFunction(x_trans) * width) +
               logfunc(roots[i], x_trans);

      if (max_psi < psi[i]) {
        max_psi = psi[i];
      }
    }

    double inner_sum = 0.0;
    for (size_t j = 0; j < psi.size(); ++j) {
      inner_sum += std::exp(psi[j] - max_psi);
    }
    logIntegral = max_psi + std::log(inner_sum);
    // --------------------------------------------------------------

    return logIntegral;
  }

  template <typename Func>
  double evaluate_sumexp_iteratively(Func const &logfunc, SingleFunction weightFunction,
                                     double a = 0.0, double b = 1.0, size_t numGaussPoints = 0,
                                     size_t incrementQuadraturePoints = 1, double tol = 1e-14) {
    auto &quadRule = base::GaussLegendreQuadRule1D::getInstance();

    double width = b - a;
    double logIntegral = 0.0;

    // performing Gauss-Legendre integration
    double err = 1e14;
    size_t iteration = 0;

    // make sure that the quadrature is done at least once
    numGaussPoints = std::min(quadRule.getMaxSupportedLevel(), numGaussPoints);

    sgpp::base::DataVector psi(numGaussPoints);
    double logIntegral_old = 0.0;
    while (err > tol && numGaussPoints <= quadRule.getMaxSupportedLevel()) {
      // --------------------------------------------------------------
      // do quadrature
      quadRule.getLevelPointsAndWeightsNormalized(numGaussPoints, roots, weights);
      psi.resize(numGaussPoints);

      double max_psi = -std::numeric_limits<double>::max();

      for (size_t i = 0; i < roots.getSize(); ++i) {
        double x_trans = a + width * roots[i];
        psi[i] = std::log(weights[i]) + std::log(weightFunction(x_trans) * width) +
                 logfunc(roots[i], x_trans);

        if (max_psi < psi[i]) {
          max_psi = psi[i];
        }
      }

      double inner_sum = 0.0;
      for (size_t j = 0; j < psi.size(); ++j) {
        inner_sum += std::exp(psi[j] - max_psi);
      }
      logIntegral = max_psi + std::log(inner_sum);
      // --------------------------------------------------------------
      if (iteration > 0) {
        err = std::fabs(logIntegral_old - logIntegral);
      }

      logIntegral_old = logIntegral;
      numGaussPoints += incrementQuadraturePoints;
      iteration += 1;

      if (incrementQuadraturePoints == 0) {
        break;
      }
    }

    return logIntegral;
  }

  template <typename Func>
  double evaluate_logsumexp_iteratively(Func const &logfunc, SingleFunction weightFunction,
                                        double a = 0.0, double b = 1.0, size_t numGaussPoints = 0,
                                        size_t incrementQuadraturePoints = 1, double tol = 1e-14) {
    return std::log(evaluate_sumexp_iteratively(logfunc, weightFunction, a, b, numGaussPoints,
                                                incrementQuadraturePoints, tol));
  }
};

} /* namespace combigrid */
} /* namespace sgpp */

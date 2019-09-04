// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/definitions.hpp>

#include <cmath>
#include <vector>

namespace sgpp {

namespace combigrid {
// ----------------------------------------------------------------------------------
struct Ishigami {
  double tolerance = 5e-4;
  size_t numDims = 3;
  double pi_4 = M_PI * M_PI * M_PI * M_PI;
  static constexpr double a = 7.0;
  static constexpr double b = 0.1;
  double variance = a * a / 8. + b * pi_4 / 5 + b * b * pi_4 * pi_4 / 18. + 0.5;
  double mean = 3.5;
  std::vector<double> sobolIndices{0.3138, 0.4424, 0.0, 0.0, 0.2436, 0.0, 0.0};
  std::vector<double> totalSobolIndices{0.5574, 0.4424, 0.2436};
  std::vector<double> bounds{0, 1};

  static double eval(sgpp::base::DataVector const& v) {
    // transform [0, 1] -> [-pi, pi]
    sgpp::base::DataVector x(v);
    x.mult(2 * M_PI);
    sgpp::base::DataVector pi_dvec(v.getSize(), M_PI);
    x.sub(pi_dvec);

    // evaluate the Ishigami function
    return std::sin(x[0]) + a * std::sin(x[1]) * std::sin(x[1]) +
           b * x[2] * x[2] * x[2] * x[2] * std::sin(x[0]);
  }
};

// ----------------------------------------------------------------------------------
struct Parabola {
  double tolerance = 1e-13;
  size_t numDims = 2;

  double alpha1 = 5.0;
  double beta1 = 4.0;

  double alpha2 = 3.0;
  double beta2 = 2.0;

  double c1 = std::tgamma(alpha1 + beta1) / (std::tgamma(alpha1) * std::tgamma(beta1));
  double c2 = std::tgamma(alpha2 + beta2) / (std::tgamma(alpha2) * std::tgamma(beta2));

  double mean = c1 * c2 / 4725.0;
  double variance = c1 * c1 * c1 * c2 * c2 * c2 / 75014100000.0 -
                    2 * c1 * c1 * c2 * c2 / 22325625.0 + 4.0 * c1 * c2 / 24255.0;
  std::vector<double> bounds{0, 1};

  static double eval(sgpp::base::DataVector const& v) {
    double ans = 1.0;
    for (size_t idim = 0; idim < v.getSize(); idim++) {
      ans *= 4.0 * v[idim] * (1.0 - v[idim]);
    }
    return ans;
  }
};
// ----------------------------------------------------------------------------------
struct Parabola_uniform {
  static double mean(size_t numDims) { return std::pow(2. / 3., numDims); }

  static double variance(size_t numDims) {
    return std::pow(16. / 30., numDims) - std::pow(mean(numDims), 2);
  }

  static void bounds(size_t numDims, std::vector<double>& bounds) {
    bounds.resize(2 * numDims);
    for (size_t i = 0; i < bounds.size(); i++) {
      if (i % 2 == 0) {
        bounds[i] = 0.0;
      } else {
        bounds[i] = 1.0;
      }
    }
  }

  static double eval(sgpp::base::DataVector const& x) {
    double ans = 1.0;
    for (size_t idim = 0; idim < x.getSize(); idim++) {
      ans *= 4 * x[idim] * (1.0 - x[idim]);
    }
    return ans;
  }
};
// ----------------------------------------------------------------------------------
struct CO2 {
  double tolerance = 1e-12;
  size_t numDims = 1;

  double logmean = std::log(1e-12);
  double stddev = std::exp(-1);
  double mean = 0.894066628227;
  double variance = 0.0102415510639;
  std::vector<double> bounds{3.87672392696e-13, 2.57949758311e-12};

  static double eval(sgpp::base::DataVector const& v) {
    return std::sin(v[0] * v[0]) + std::cos(2 * v[0]);
  }
};
// ----------------------------------------------------------------------------------
struct AtanUniform {
  size_t numDims = 2;

  double mean = 3.514491266446367;
  double variance = 3.453103593932221;

  static double eval(sgpp::base::DataVector const& v) {
    return std::atan(50.0 * (v[0] - 0.35)) + M_PI / 2.0 + 4.0 * std::pow(v[1], 3.0) +
           std::exp(v[0] * v[1] - 1.0);
  }
};
// ----------------------------------------------------------------------------------
struct AtanBeta {
  size_t numDims = 2;

  double alpha1 = 5.0;
  double beta1 = 4.0;

  double alpha2 = 3.0;
  double beta2 = 2.0;

  std::vector<double> bounds{0, 1};

  double mean = 0.0;
  double variance = 0.0;

  static double eval(sgpp::base::DataVector const& v) {
    return std::atan(50.0 * (v[0] - 0.35)) + M_PI / 2.0 + 4.0 * std::pow(v[1], 3.0) +
           std::exp(v[0] * v[1] - 1.0);
  }
};
// ----------------------------------------------------------------------------------

struct Debugfct {
  size_t numDims = 1;
  std::vector<double> bounds{0, 1};
  double mean = 0.5;
  double variance = 0.083333333333333333333;
  static double eval(sgpp::base::DataVector const& v) { return v[0]; }
};
// ----------------------------------------------------------------------------------
struct Genz {
  static void bounds(size_t numDims, std::vector<double>& bounds) {
    bounds.resize(2 * numDims);
    for (size_t i = 0; i < bounds.size(); i++) {
      if (i % 2 == 0) {
        bounds[i] = 0.0;
      } else {
        bounds[i] = 1.0;
      }
    }
  }

  static constexpr double w = 0.0;

  static double eval(sgpp::base::DataVector const& x) {
    double ans = 2.0 * M_PI * w;
    size_t numDims = x.size();
    for (size_t k = 0; k < numDims; k++) {
      ans += 4.5 * (static_cast<double>(k) + 0.5) / static_cast<double>(numDims) * x[k];
    }
    return std::cos(ans);
  }
};
// ----------------------------------------------------------------------------------

}  // namespace combigrid
}  // namespace sgpp

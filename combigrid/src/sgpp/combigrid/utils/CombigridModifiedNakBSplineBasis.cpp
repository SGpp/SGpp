// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "CombigridModifiedNakBSplineBasis.hpp"
#include <algorithm>
#include <vector>
#include "CombigridBSplineBasis.hpp"

namespace sgpp {
namespace combigrid {

// log2 for size_t input and output
size_t logarithm2(size_t n) { return ((n <= 2) ? 1 : 1 + logarithm2(n / 2)); }

/**
 * Modified Not-a-Knot B spline basis
 */
double expUniformModifiedNakBspline(double const& x, size_t const& degree, size_t i,
                                    std::vector<double> points) {
  // derive level from number of gridpoints
  size_t l = combigrid::logarithm2(points.size() + 1);

  const size_t hInv = 1 << l;
  double t = x * static_cast<double>(hInv) - static_cast<double>(i);

  //  std::cout << l << " " << i << " " << hInv << " " << x << std::endl;
  switch (degree) {
    case 3:
      if (l == 1) {
        // l = 1, i = 0
        return 1.0;
      } else if ((i > 0) && (i < hInv - 2)) {
        // l >= 3, 0 < i < 2^l - 2	(the unmodified cases i=0 and i=2^l-1)
        points.insert(points.begin(), 1, 0);  // expensive!
        points.push_back(1);
        return sgpp::combigrid::expUniformNakBspline(x, degree, i + 1, points);
      } else {
        if (i > (hInv / 2 - 1)) {
          i = hInv - i - 2;
          //          t *= -1.0;
          t = (1 - x) * static_cast<double>(hInv) - static_cast<double>(i);
        }

        if (l == 2) {
          // l = 2, i = 0
          if ((t < -1.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 1.0;
            double result = 1.0 / 40.0;
            result *= t;
            result = -3.0 / 5.0 + result * t;
            result = 6.0 / 5.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 40.0;
            result = 3.0 / 20.0 + result * t;
            result = -3.0 / 10.0 + result * t;
            result = 1.0 / 5.0 + result * t;
            return result;
          }
        } else {
          // l >= 3, i = 0
          if ((t < -1.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 1.0;
            double result = 1.0 / 24.0;
            result *= t;
            result = -3.0 / 4.0 + result * t;
            result = 5.0 / 4.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 12.0;
            result = 1.0 / 4.0 + result * t;
            result = -1.0 / 4.0 + result * t;
            result = 1.0 / 12.0 + result * t;
            return result;
          }
        }
      }

      // cases degree = 1, 5 and 7 could be copied from NakBsplineModifiedBasis and adjusted
      // Also there are no quadrature routines for these basis functions so they can only be used
      // for interpolation
    default:
      return 0.0;
  }
}
}  // namespace combigrid
}  // namespace sgpp

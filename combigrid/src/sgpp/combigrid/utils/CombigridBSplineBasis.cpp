// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/utils/CombigridBSplineBasis.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

size_t log2(size_t n) { return ((n < 2) ? 1 : 1 + static_cast<size_t>(log2(n / 2))); }

/**
 * Nodal Not-a-Knot B spline basis
 * always constant on level 0
 * degree 1 : hats from level 1 onwards
 * degree 3 : polynomials on level 1, Not-a-Knot B-splines from level 2 onwards
 * degree 5 : polynomials on level 1,2, Not-a-KNot B-Splines from level 3 onwards
 *
 * the level is determined by the length of the points vector containing the grid points
 */
double expUniformNakBspline(double const& x, size_t const& degree, size_t i,
                            std::vector<double> const& points) {
  // derive level from number of gridpoints
  size_t l = 0;
  if (points.size() > 1) {
    l = log2(points.size() - 2);
  }
  const size_t hInv = 1 << l;
  double t = x * static_cast<double>(hInv) - static_cast<double>(i);

  switch (degree) {
    case 1:
      if (l == 0) {
        return 1.0;
      }
      return std::max(1.0 - std::abs(t), 0.0);

    // degree 3: polynomials on Level 0 and 1, nak Bsplines from Level 2 on
    case 3: {
      sgpp::base::BsplineBasis<size_t, size_t> bsplineBasis(3);
      if (l == 0) {
        // l = 0, i = 0
        return 1.0;

      } else if (l == 1) {
        // Lagrange polynomials
        if (i == 0) {
          // l = 1, i = 0
          return (2 * x - 3) * x + 1;
        } else if (i == 1) {
          // l = 1, i = 1
          return 4 * x * (1 - x);
        } else {
          // l = 1, i = 2
          return (2 * x - 1) * x;
        }
      } else if ((i > 3) && (i < hInv - 3)) {
        // l >= 4, 3 < i < 2^l - 3
        return bsplineBasis.eval(l, i, x);
      } else {
        if (i > hInv / 2) {
          i = hInv - i;
          t *= -1.0;
        }

        if (i == 0) {
          // l >= 2, i = 0
          if ((t < 0.0) || (t > 2.0)) {
            return 0.0;
          } else {
            double result = -4.1666666666666664e-02;
            result = 2.5000000000000000e-01 + result * t;
            result = -5.0000000000000000e-01 + result * t;
            result = 3.3333333333333331e-01 + result * t;
            return result;
          }
        } else if ((l == 2) && (i == 1)) {
          // l = 2, i = 1
          if ((t < -1.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 1.0;
            double result = 1.0 / 10.0;
            result = -9.0 / 20.0 + result * t;
            result = 3.0 / 10.0 + result * t;
            result = 3.0 / 5.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 40.0;
            result = 3.0 / 20.0 + result * t;
            result = -3.0 / 10.0 + result * t;
            result = 1.0 / 5.0 + result * t;
            return result;
          }

        } else if (l == 2) {
          // l = 2, i = 2
          if ((t < -2.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 2.0;
            double result = -8.3333333333333329e-02;
            result = 2.0000000000000001e-01 + result * t;
            result = 2.0000000000000001e-01 + result * t;
            result = 6.6666666666666666e-02 + result * t;
            return result;
          } else {
            double result = 8.3333333333333329e-02;
            result = -2.9999999999999999e-01 + result * t;
            result *= t;
            result = 5.9999999999999998e-01 + result * t;
            return result;
          }
        } else if (i == 1) {
          // l >= 3, i = 1
          if ((t < -1.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 1.0;
            double result = 1.0 / 8.0;
            result = -1.0 / 2.0 + result * t;
            result = 1.0 / 4.0 + result * t;
            result = 7.0 / 12.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 12.0;
            result = 1.0 / 4.0 + result * t;
            result = -1.0 / 4.0 + result * t;
            result = 1.0 / 12.0 + result * t;
            return result;
          }
        } else if (i == 2) {
          // l >= 3, i = 2
          if ((t < -2.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 2.0;
            double result = -1.2500000000000000e-01;
            result = 2.5000000000000000e-01 + result * t;
            result = 2.5000000000000000e-01 + result * t;
            result = 8.3333333333333329e-02 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = 2.9166666666666669e-01;
            result = -5.0000000000000000e-01 + result * t;
            result = -2.5000000000000000e-01 + result * t;
            result = 5.8333333333333337e-01 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.2500000000000000e-01;
            result = 3.7500000000000000e-01 + result * t;
            result = -3.7500000000000000e-01 + result * t;
            result = 1.2500000000000000e-01 + result * t;
            return result;
          }
        } else {
          // l >= 3, i = 3
          if ((t < -3.0) || (t > 2.0)) {
            return 0.0;
          } else if (t < -1.0) {
            t += 3.0;
            double result = 1.0 / 24.0;
            result *= t;
            result *= t;
            result *= t;
            return result;
          } else if (t < 0.0) {
            t += 1.0;
            double result = -3.0 / 8.0;
            result = 1.0 / 4.0 + result * t;
            result = 1.0 / 2.0 + result * t;
            result = 1.0 / 3.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = 11.0 / 24.0;
            result = -7.0 / 8.0 + result * t;
            result = -1.0 / 8.0 + result * t;
            result = 17.0 / 24.0 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = -1.0 / 6.0;
            result = 1.0 / 2.0 + result * t;
            result = -1.0 / 2.0 + result * t;
            result = 1.0 / 6.0 + result * t;
            return result;
          }
        }
      }
    }
    // degree 5: Levels 0,1 and 2 polynomials, nak B splines from Level 3 on
    case 5: {
      sgpp::base::BsplineBasis<size_t, size_t> bsplineBasis(5);
      if (l == 0) {
        // l = 0, i = 0
        return 1.0;

      } else if (l == 1) {
        // Lagrange polynomials
        if (i == 0) {
          // l = 1, i = 0
          return (2 * x - 3) * x + 1;
        } else if (i == 1) {
          // l = 1, i = 1
          return 4 * x * (1 - x);
        } else {
          // l = 1, i = 2
          return (2 * x - 1) * x;
        }
      } else if (l == 2) {
        if (i == 0) {
          // l = 2, i = 0
          return 1.0 / 3.0 * (x - 1) * (2 * x - 1) * (4 * x - 3) * (4 * x - 1);
        } else if (i == 1) {
          // l = 2, i = 1
          return -16.0 / 3.0 * x * (2 * x - 1) * (4 * x - 3) * (x - 1);
        } else if (i == 2) {
          // l = 2, i = 2
          return 4 * x * (4 * x - 1) * (4 * x - 3) * (x - 1);
        } else if (i == 3) {
          // l = 2, i = 3
          return -16.0 / 3.0 * (x - 1) * x * (2 * x - 1) * (4 * x - 1);
        } else {
          // l = 2, i = 4
          return 1.0 / 3.0 * x * (4 * x - 1) * (2 * x - 1) * (4 * x - 3);
        }
      } else if ((i > 5) && (i < hInv - 5)) {
        // l >= 4, 5 < i < 2^l - 5
        return bsplineBasis.eval(l, i, x);
      } else {
        if (i > hInv / 2) {
          i = hInv - i;
          t *= -1.0;
        }

        if ((l == 3) && (i == 3)) {
          // l = 3, i = 3
          if ((t < -3.0) || (t > 5.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 3.0;
            double result = 107.0 / 30240.0;
            result = -17.0 / 756.0 + result * t;
            result = 1.0 / 378.0 + result * t;
            result = 37.0 / 378.0 + result * t;
            result = 109.0 / 756.0 + result * t;
            result = 253.0 / 3780.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -397.0 / 30240.0;
            result = 185.0 / 6048.0 + result * t;
            result = 155.0 / 3024.0 + result * t;
            result = -415.0 / 3024.0 + result * t;
            result = -1165.0 / 6048.0 + result * t;
            result = 2965.0 / 6048.0 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 233.0 / 30240.0;
            result = -53.0 / 1512.0 + result * t;
            result = 8.0 / 189.0 + result * t;
            result = 13.0 / 189.0 + result * t;
            result = -97.0 / 378.0 + result * t;
            result = 433.0 / 1890.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 4320.0;
            result = 1.0 / 288.0 + result * t;
            result = -1.0 / 48.0 + result * t;
            result = 1.0 / 16.0 + result * t;
            result = -3.0 / 32.0 + result * t;
            result = 9.0 / 160.0 + result * t;
            return result;
          }
        } else if ((l == 3) && (i == 4)) {
          // l=3, i=4
          if ((t < -4.0) || (t > 4.0)) {
            return 0.0;
          } else if (t < -1.0) {
            t += 4.0;
            double result = -1.3888888888888889e-03;
            result = 4.6296296296296294e-03 + result * t;
            result = 9.2592592592592587e-03 + result * t;
            result = 9.2592592592592587e-03 + result * t;
            result = 4.6296296296296294e-03 + result * t;
            result = 9.2592592592592596e-04 + result * t;
            return result;
          } else if (t < 0.0) {
            t += 1.0;
            double result = 1.2500000000000001e-02;
            result = -1.6203703703703703e-02 + result * t;
            result = -6.0185185185185182e-02 + result * t;
            result = -3.2407407407407406e-02 + result * t;
            result = 2.4768518518518517e-01 + result * t;
            result = 3.8564814814814813e-01 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -1.2500000000000001e-02;
            result = 4.6296296296296294e-02 + result * t;
            result *= t;
            result = -1.8518518518518517e-01 + result * t;
            result *= t;
            result = 5.3703703703703709e-01 + result * t;
            return result;
          } else {
            t -= 1.0;
            double result = 1.3888888888888889e-03;
            result = -1.6203703703703703e-02 + result * t;
            result = 6.0185185185185182e-02 + result * t;
            result = -3.2407407407407406e-02 + result * t;
            result = -2.4768518518518517e-01 + result * t;
            result = 3.8564814814814813e-01 + result * t;
            return result;
          }
        } else if (i == 0) {
          // l>=3, i=0
          if ((t < 0.0) || (t > 3.0)) {
            return 0.0;
          } else {
            double result = -3.9682539682539683e-04;
            result = 5.9523809523809521e-03 + result * t;
            result = -3.5714285714285712e-02 + result * t;
            result = 1.0714285714285714e-01 + result * t;
            result = -1.6071428571428573e-01 + result * t;
            result = 9.6428571428571433e-02 + result * t;
            return result;
          }
        } else if (i == 1) {
          // l >= 3, i = 1
          if ((t < -1.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 2.0) {
            t += 1.0;
            double result = 1.0 / 504.0;
            result = -1.0 / 42.0 + result * t;
            result = 2.0 / 21.0 + result * t;
            result = -2.0 / 21.0 + result * t;
            result = -5.0 / 21.0 + result * t;
            result = 47.0 / 105.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 840.0;
            result = 1.0 / 168.0 + result * t;
            result = -1.0 / 84.0 + result * t;
            result = 1.0 / 84.0 + result * t;
            result = -1.0 / 168.0 + result * t;
            result = 1.0 / 840.0 + result * t;
            return result;
          }
        } else if (i == 2) {
          // l>=3, i=2
          if ((t < -2.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 1.0) {
            t += 2.0;
            double result = -3.9682539682539680e-03;
            result = 3.5714285714285712e-02 + result * t;
            result = -7.1428571428571425e-02 + result * t;
            result = -1.1904761904761904e-01 + result * t;
            result = 2.5000000000000000e-01 + result * t;
            result = 3.8809523809523810e-01 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 7.1428571428571426e-03;
            result = -2.3809523809523808e-02 + result * t;
            result *= t;
            result = 9.5238095238095233e-02 + result * t;
            result = -1.4285714285714285e-01 + result * t;
            result = 6.6666666666666666e-02 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -2.3809523809523812e-03;
            result = 1.1904761904761904e-02 + result * t;
            result = -2.3809523809523808e-02 + result * t;
            result = 2.3809523809523808e-02 + result * t;
            result = -1.1904761904761904e-02 + result * t;
            result = 2.3809523809523812e-03 + result * t;
            return result;
          }
        } else if (i == 3) {
          // l >= 4, i = 3
          if ((t < -3.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < 0.0) {
            t += 3.0;
            double result = 1.0 / 252.0;
            result = -1.0 / 42.0 + result * t;
            result *= t;
            result = 2.0 / 21.0 + result * t;
            result = 1.0 / 7.0 + result * t;
            result = 1.0 / 15.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -23.0 / 1260.0;
            result = 1.0 / 28.0 + result * t;
            result = 1.0 / 14.0 + result * t;
            result = -5.0 / 42.0 + result * t;
            result = -1.0 / 4.0 + result * t;
            result = 163.0 / 420.0 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 19.0 / 1260.0;
            result = -1.0 / 18.0 + result * t;
            result = 2.0 / 63.0 + result * t;
            result = 8.0 / 63.0 + result * t;
            result = -2.0 / 9.0 + result * t;
            result = 34.0 / 315.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 252.0;
            result = 5.0 / 252.0 + result * t;
            result = -5.0 / 126.0 + result * t;
            result = 5.0 / 126.0 + result * t;
            result = -5.0 / 252.0 + result * t;
            result = 1.0 / 252.0 + result * t;
            return result;
          }
        } else if (i == 4) {
          // l>=4 i=4
          if ((t < -4.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < -1.0) {
            t += 4.0;
            double result = -1.9841269841269840e-03;
            result = 5.9523809523809521e-03 + result * t;
            result = 1.1904761904761904e-02 + result * t;
            result = 1.1904761904761904e-02 + result * t;
            result = 5.9523809523809521e-03 + result * t;
            result = 1.1904761904761906e-03 + result * t;
            return result;
          } else if (t < 0.0) {
            t += 1.0;
            double result = 2.5793650793650792e-02;
            result = -2.3809523809523808e-02 + result * t;
            result = -9.5238095238095233e-02 + result * t;
            result = -9.5238095238095233e-02 + result * t;
            result = 2.3809523809523808e-01 + result * t;
            result = 4.4761904761904764e-01 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -4.0873015873015874e-02;
            result = 1.0515873015873016e-01 + result * t;
            result = 6.7460317460317457e-02 + result * t;
            result = -2.6587301587301587e-01 + result * t;
            result = -2.0436507936507936e-01 + result * t;
            result = 4.9722222222222223e-01 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 2.5793650793650792e-02;
            result = -9.9206349206349201e-02 + result * t;
            result = 7.9365079365079361e-02 + result * t;
            result = 1.5873015873015872e-01 + result * t;
            result = -3.1746031746031744e-01 + result * t;
            result = 1.5873015873015872e-01 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -5.9523809523809521e-03;
            result = 2.9761904761904760e-02 + result * t;
            result = -5.9523809523809521e-02 + result * t;
            result = 5.9523809523809521e-02 + result * t;
            result = -2.9761904761904760e-02 + result * t;
            result = 5.9523809523809521e-03 + result * t;
            return result;
          }
        } else {
          // l >= 4, i = 5
          if ((t < -5.0) || (t > 3.0)) {
            return 0.0;
          } else if (t < -2.0) {
            t += 5.0;
            double result = 1.0 / 2520.0;
            result *= t;
            result *= t;
            result *= t;
            result *= t;
            result *= t;
            return result;
          } else if (t < -1.0) {
            t += 2.0;
            double result = -11.0 / 504.0;
            result = 1.0 / 168.0 + result * t;
            result = 1.0 / 28.0 + result * t;
            result = 3.0 / 28.0 + result * t;
            result = 9.0 / 56.0 + result * t;
            result = 27.0 / 280.0 + result * t;
            return result;
          } else if (t < 0.0) {
            t += 1.0;
            double result = 31.0 / 504.0;
            result = -13.0 / 126.0 + result * t;
            result = -10.0 / 63.0 + result * t;
            result = 2.0 / 63.0 + result * t;
            result = 25.0 / 63.0 + result * t;
            result = 121.0 / 315.0 + result * t;
            return result;
          } else if (t < 1.0) {
            double result = -181.0 / 2520.0;
            result = 103.0 / 504.0 + result * t;
            result = 11.0 / 252.0 + result * t;
            result = -113.0 / 252.0 + result * t;
            result = -61.0 / 504.0 + result * t;
            result = 1543.0 / 2520.0 + result * t;
            return result;
          } else if (t < 2.0) {
            t -= 1.0;
            double result = 11.0 / 280.0;
            result = -13.0 / 84.0 + result * t;
            result = 1.0 / 7.0 + result * t;
            result = 4.0 / 21.0 + result * t;
            result = -3.0 / 7.0 + result * t;
            result = 23.0 / 105.0 + result * t;
            return result;
          } else {
            t -= 2.0;
            double result = -1.0 / 120.0;
            result = 1.0 / 24.0 + result * t;
            result = -1.0 / 12.0 + result * t;
            result = 1.0 / 12.0 + result * t;
            result = -1.0 / 24.0 + result * t;
            result = 1.0 / 120.0 + result * t;
            return result;
          }
        }
      }
    }

    default:
      return 0.0;
  }
}

double nonUniformBSpline(double const& x, size_t const& deg, size_t const& index,
                         std::vector<double> const& xi) {
  if (deg == 0) {
    // characteristic function of [xi[k], xi[k+1])
    return (((x >= xi[index]) && (x < xi[index + 1])) ? 1.0 : 0.0);
  } else if ((x < xi[index]) || (x >= xi[index + deg + 1])) {
    // out of support
    return 0.0;
  } else {
    // Cox-de-Boor recursion
    return (x - xi[index]) / (xi[index + deg] - xi[index]) *
               nonUniformBSpline(x, deg - 1, index, xi) +
           (1.0 - (x - xi[index + 1]) / (xi[index + deg + 1] - xi[index + 1])) *
               nonUniformBSpline(x, deg - 1, index + 1, xi);
  }
}

double LagrangePolynomial(double const& x, std::vector<double> const& xValues, size_t const& k) {
  double res = 1.0;
  for (size_t m = 0; m < xValues.size(); m++) {
    if (k != m) {
      res *= (x - xValues[m]) / (xValues[k] - xValues[m]);
    }
  }
  return res;
}
}  // namespace combigrid
}  // namespace sgpp

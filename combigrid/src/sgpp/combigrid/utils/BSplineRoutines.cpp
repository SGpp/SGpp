// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/utils/BSplineRoutines.hpp>
#include <vector>

constexpr size_t log2(size_t n) { return ((n < 2) ? 1 : 1 + log2(n / 2)); }

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

    // degree 3: polynomials on Level 0 and 1, nak Bsplines from Level 3 on
    case 3: {
      sgpp::base::BsplineBasis<size_t, size_t> bsplineBasis(3);
      if (l == 0) {
        if (i == 0) {
          // l = 0, i = 0
          return 1.0;
        }

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
        if (i == 0) {
          // l = 0, i = 0
          return 1.0;
        }

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
        } else if (i == 4) {
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

// ToDo (rehmemk) Test other strategies for outer points for example placing them uniform using
// the
// distance of the last point inside to the boundary

std::vector<double> createdeg1Knots(std::vector<double> const& xValues) {
  size_t degree = 1;

  // this offset works only for odd degrees. if even degrees shall be supported it must be
  // generalized
  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * offset, 0);
  xi.insert(xi.begin() + offset, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset; i++) {
    xi[offset - i - 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }
  return xi;
}

std::vector<double> createdeg3NakKnots(std::vector<double> const& xValues) {
  size_t degree = 3;
  // On levels 0 and 1 only Lagrange polynomials and not nak B splines are used. no need for extra
  // points
  if (xValues.size() < 5) {
    return xValues;
  }
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * offset + 2, 0);

  xi.insert(xi.begin() + offset + 1, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 1; i++) {
    xi[offset - i] = -xValues.at(i + 1) + 2 * xValues[0];
    xi[xValues.size() + offset + i + 1] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues.at(xValues.size() - i - 2));
  }

  xi.erase(xi.begin() + offset + 2);
  xi.erase(xi.end() - offset - 3);
  return xi;
}

std::vector<double> createdeg5NakKnots(std::vector<double> const& xValues) {
  size_t degree = 5;
  // On levels 0,1 and 2 only Lagrange polynomials and not nak B splines are used. no need for
  // extra
  // points
  if (xValues.size() < 9) {
    return xValues;
  }
  // create a vector xi that holds the gridpoints and continues to the left and right by mirroring
  // at 0 and 1

  size_t offset = (degree + 1) / 2;
  std::vector<double> xi(2 * (offset + 2), 0);

  xi.insert(xi.begin() + offset + 2, xValues.begin(), xValues.end());

  for (size_t i = 0; i < offset + 2; i++) {
    xi[offset - i + 1] = -xValues[i + 1] + 2 * xValues[0];
    xi[xValues.size() + offset + i + 2] =
        xValues[xValues.size() - 1] +
        (xValues[xValues.size() - 1] - xValues[xValues.size() - i - 2]);
  }

  xi.erase(xi.begin() + offset + 3);
  xi.erase(xi.begin() + offset + 3);
  xi.erase(xi.end() - offset - 4);
  xi.erase(xi.end() - offset - 4);
  return xi;
}

std::vector<double> createNakKnots(std::vector<double> const& xValues, size_t const& degree) {
  if (degree == 1) {
    return createdeg1Knots(xValues);
  } else if (degree == 3) {
    return createdeg3NakKnots(xValues);
  } else if (degree == 5) {
    return createdeg5NakKnots(xValues);
  } else {
    throw std::invalid_argument("BSplineRoutines: only B-Spline degrees 1,3 and 5 supported");
  }
}

// ToDo (rehmemk) use unidirectional principle instead of global SLE solving to speed this up!

sgpp::combigrid::GridFunction BSplineCoefficientGridFunction(
    sgpp::combigrid::MultiFunction func, sgpp::combigrid::CombiHierarchies::Collection grids,
    size_t degree) {
  sgpp::combigrid::GridFunction gf([func, grids,
                                    degree](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
    size_t numDimensions = grid->getDimension();
    // stores the values of the objective function
    auto funcStorage = std::make_shared<sgpp::combigrid::CombigridTreeStorage>(grids, func);

    // ToDo (rehmemk) B-spline interpolation and B-spline quadrature can be mixed (one in one
    // dimension the other in another dimension and so on). Combining B-splines and other basis
    // functions has not been tested yet.

    sgpp::combigrid::CombiEvaluators::Collection interpolEvaluators(
        numDimensions, sgpp::combigrid::CombiEvaluators::BSplineInterpolation(degree));

    auto coefficientTree = std::make_shared<sgpp::combigrid::TreeStorage<double>>(numDimensions);
    auto level = grid->getLevel();
    std::vector<size_t> numGridPointsVec = grid->numPoints();
    size_t numGridPoints = 1;
    for (size_t i = 0; i < numGridPointsVec.size(); i++) {
      numGridPoints *= numGridPointsVec[i];
    }

    std::vector<bool> orderingConfiguration;

    sgpp::combigrid::CombiEvaluators::Collection evalCopy(numDimensions);
    for (size_t dim = 0; dim < numDimensions; ++dim) {
      evalCopy[dim] = interpolEvaluators[dim]->cloneLinear();
      bool needsSorted = evalCopy[dim]->needsOrderedPoints();
      auto gridPoints = grids[dim]->getPoints(level[dim], needsSorted);
      orderingConfiguration.push_back(needsSorted);
      evalCopy[dim]->setGridPoints(gridPoints);
    }
    sgpp::base::DataMatrix A(numGridPoints, numGridPoints);
    sgpp::base::DataVector coefficients_sle(numGridPoints);
    sgpp::base::DataVector functionValues(numGridPoints);

    // Creates an iterator that yields the multi-indices of all grid points in the grid.
    sgpp::combigrid::MultiIndexIterator it(grid->numPoints());

    auto funcIter = funcStorage->getGuidedIterator(level, it, orderingConfiguration);

    for (size_t ixEvalPoints = 0; funcIter->isValid(); ++ixEvalPoints, funcIter->moveToNext()) {
      auto gridPoint = grid->getGridPoint(funcIter->getMultiIndex());
      functionValues[ixEvalPoints] = funcIter->value();

      std::vector<std::vector<double>> basisValues;
      for (size_t dim = 0; dim < numDimensions; ++dim) {
        evalCopy[dim]->setParameter(sgpp::combigrid::FloatScalarVector(gridPoint[dim]));
        auto basisValues1D = evalCopy[dim]->getBasisValues();
        // basis values at gridPoint
        std::vector<double> basisValues1D_vec(basisValues1D.size());
        for (size_t i = 0; i < basisValues1D.size(); i++) {
          basisValues1D_vec[i] = basisValues1D[i].value();
        }
        basisValues.push_back(basisValues1D_vec);
      }

      sgpp::combigrid::MultiIndexIterator innerIter(grid->numPoints());
      for (size_t ixBasisFunctions = 0; innerIter.isValid();
           ++ixBasisFunctions, innerIter.moveToNext()) {
        double splineValue = 1.0;
        auto innerIndex = innerIter.getMultiIndex();
        for (size_t dim = 0; dim < numDimensions; ++dim) {
          splineValue *= basisValues[dim][innerIndex[dim]];
        }
        A.set(ixEvalPoints, ixBasisFunctions, splineValue);
      }
    }

    sgpp::optimization::FullSLE sle(A);

    //    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
    //    bool solved = 0;
    //    sgpp::combigrid::Stopwatch watch;
    //    watch.start();
    //    if (numGridPoints < 500) {
    //      sgpp::optimization::sle_solver::Armadillo solver;
    //      solved = solver.solve(sle, functionValues, coefficients_sle);
    //    } else {
    //      sgpp::optimization::sle_solver::UMFPACK solver;
    //      solved = solver.solve(sle, functionValues, coefficients_sle);
    //    }
    //    std::cout << numGridPoints << " " << watch.elapsedSeconds() << std::endl;

    sgpp::optimization::sle_solver::UMFPACK solver;
    sgpp::optimization::Printer::getInstance().setVerbosity(-1);
    bool solved = solver.solve(sle, functionValues, coefficients_sle);

    /*
    std::cout << A.toString() << std::endl;
    std::cout << "fct: ";
    for (size_t i = 0; i < functionValues.size(); i++) {
      std::cout << functionValues[i] << " ";
    }
    std::cout << "\ncoeff: ";
    for (size_t i = 0; i < coefficients_sle.size(); i++) {
      std::cout << coefficients_sle[i] << " ";
    }
    std::cout << "\n";
    std::cout << "--------" << std::endl;
    */

    if (!solved) {
      exit(-1);
    }

    it.reset();
    for (size_t vecIndex = 0; it.isValid(); ++vecIndex, it.moveToNext()) {
      coefficientTree->set(it.getMultiIndex(), coefficients_sle[vecIndex]);
    }

    return coefficientTree;
  });
  return gf;
}

// levelManager must be an AveragingLevelManager. Otherwise this makes no sense
std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> createBsplineVarianceRefinementOperation(
    size_t degree, size_t numDimensions, sgpp::combigrid::MultiFunction func,
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);

  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));

  //  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
  //      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;

  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);
  return Operation;
}

std::shared_ptr<sgpp::combigrid::CombigridMultiOperation>
createBsplineVarianceRefinementOperationWithWeightsAndBounds(
    size_t degree, sgpp::combigrid::MultiFunction func,
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager,
    sgpp::combigrid::WeightFunctionsCollection weightFunctionsCollection,
    sgpp::base::DataVector bounds) {
  size_t numDimensions = weightFunctionsCollection.size();
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineScalarProduct, degree);

  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));

  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::VARIANCE;

  for (size_t d = 0; d < numDimensions; d++) {
    Evaluators[d]->setWeightFunction(weightFunctionsCollection[d]);
    Evaluators[d]->setBounds(bounds[2 * d], bounds[2 * d + 1]);
  }

  bool exploitNesting = false;
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);

  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);

  return Operation;
}

std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> createBsplineLinearRefinementOperation(
    size_t degree, size_t numDimensions, sgpp::combigrid::MultiFunction func,
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager) {
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::GridFunction gf = BSplineCoefficientGridFunction(func, pointHierarchies, degree);

  sgpp::combigrid::EvaluatorConfiguration EvalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, degree);

  sgpp::combigrid::CombiEvaluators::MultiCollection Evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(EvalConfig));

  //  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
  //      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType auxiliarySummationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;

  bool exploitNesting = false;
  auto Operation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, Evaluators, levelManager, gf, exploitNesting,
      auxiliarySummationStrategyType);
  return Operation;
}

std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> createBsplineLinearCoefficientOperation(
    size_t degree, size_t numDimensions,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage) {
  //    sgpp::combigrid::MultiFunction func) {
  sgpp::combigrid::EvaluatorConfiguration evalConfig(
      sgpp::combigrid::CombiEvaluatorTypes::Multi_BSplineInterpolation, degree);
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  sgpp::combigrid::CombiEvaluators::MultiCollection evaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::createCombiMultiEvaluator(evalConfig));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  std::shared_ptr<sgpp::combigrid::LevelManager> dummyLevelManager(
      new sgpp::combigrid::AveragingLevelManager());
  auto interpolationOperation = std::make_shared<sgpp::combigrid::CombigridMultiOperation>(
      pointHierarchies, evaluators, dummyLevelManager, coefficientStorage, summationStrategyType);
  return interpolationOperation;
}

std::shared_ptr<sgpp::combigrid::CombigridOperation>
createexpUniformBsplineQuadratureCoefficientOperation(
    size_t degree, size_t numDimensions,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage) {
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree));
  sgpp::combigrid::CombiHierarchies::Collection pointHierarchies(
      numDimensions, sgpp::combigrid::CombiHierarchies::expUniformBoundary());
  std::shared_ptr<sgpp::combigrid::LevelManager> dummyLevelManager(
      new sgpp::combigrid::AveragingLevelManager());
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  auto quadOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, quadEvaluators, dummyLevelManager, coefficientStorage,
      summationStrategyType);
  return quadOperation;
}

std::shared_ptr<sgpp::combigrid::CombigridOperation> createBsplineQuadratureCoefficientOperation(
    size_t degree, size_t numDimensions,
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager,
    sgpp::combigrid::CombiHierarchies::Collection pointHierarchies,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage) {
  sgpp::combigrid::CombiEvaluators::Collection quadEvaluators(
      numDimensions, sgpp::combigrid::CombiEvaluators::BSplineQuadrature(degree));
  sgpp::combigrid::FullGridSummationStrategyType summationStrategyType =
      sgpp::combigrid::FullGridSummationStrategyType::LINEAR;
  auto quadOperation = std::make_shared<sgpp::combigrid::CombigridOperation>(
      pointHierarchies, quadEvaluators, levelManager, coefficientStorage, summationStrategyType);
  return quadOperation;
}
void printLevelStructure(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelstructure) {
  auto it = levelstructure->getStoredDataIterator();
  while (it->isValid()) {
    sgpp::combigrid::MultiIndex index = it->getMultiIndex();
    for (auto& i : index) {
      std::cout << i << " ";
    }
    std::cout << "\n";
    it->moveToNext();
  }
}

sgpp::base::DataMatrix convertLevelStructureToMatrix(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelstructure, size_t numDims) {
  sgpp::base::DataMatrix levelstructureMatrix(0, numDims);
  auto it = levelstructure->getStoredDataIterator();
  while (it->isValid()) {
    sgpp::combigrid::MultiIndex index = it->getMultiIndex();
    sgpp::base::DataVector row;
    for (auto& i : index) {
      row.push_back(static_cast<double>(i));
    }
    levelstructureMatrix.appendRow(row);
    it->moveToNext();
  }
  return levelstructureMatrix;
}

sgpp::base::DataMatrix convertLevelStructureToGridPoints(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    size_t numDimensions, size_t degree) {
  sgpp::base::DataMatrix gridpointMatrix(0, numDimensions);
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDimensions, degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);
  for (size_t q = 0; q < gridStorage.getSize(); q++) {
    auto point = gridStorage.getPoint(q);
    sgpp::base::DataVector row;
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      row.push_back(point.getStandardCoordinate(d));
    }
    gridpointMatrix.appendRow(row);
  }
  return gridpointMatrix;
}

void printSGGridToFile(std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
                       size_t numDimensions, size_t degree) {
  std::string plotstr = "/home/rehmemk/SGS_Sync/Plotting/combigrid_bsplines/convertedGrid.dat";
  remove(plotstr.c_str());
  std::ofstream plotfile;
  plotfile.open(plotstr.c_str(), std::ios::app);
  plotfile << "#grid points" << std::endl;
  // convert level structure to SG
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDimensions, degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);
  for (size_t q = 0; q < gridStorage.getSize(); q++) {
    auto point = gridStorage.getPoint(q);
    for (size_t d = 0; d < gridStorage.getDimension() - 1; d++) {
      plotfile << point.getStandardCoordinate(d) << ", ";
    }
    plotfile << point.getStandardCoordinate(gridStorage.getDimension() - 1);
    plotfile << "\n";
  }
  plotfile.close();
}

std::vector<double> calculateBsplineMeanAndVariance(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    size_t numDimensions, size_t degree,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage) {
  //  size_t numthreads = 4;

  // create CT interpolation operation
  auto interpolationOperation =
      createBsplineLinearCoefficientOperation(degree, numDimensions, coefficientStorage);

  // convert level structure to SG
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDimensions, degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);

  // interpolate on SG
  sgpp::base::DataMatrix interpolParams(numDimensions, gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    interpolParams.setColumn(i, p);
  }

  // obtain function values from combigrid surrogate
  interpolationOperation->setParameters(interpolParams);
  interpolationOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector f_values = interpolationOperation->getResult();

  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::UMFPACK sleSolver;
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }

  //-----------error calculation ---------------------------------
  //  sgpp::optimization::InterpolantScalarFunction u(*grid, alpha);
  //  double diff = 0;
  //  double maxErr = 0;
  //  double L2Err = 0;
  //  size_t numMCpoints = 10000;
  //  sgpp::quadrature::NaiveSampleGenerator generator(numDimensions);
  //  sgpp::base::DataVector p(numDimensions);
  //  for (size_t i = 0; i < numMCpoints; i++) {
  //    generator.getSample(p);
  //    double functionevaluation = 1;  // f(p)
  //    diff = fabs(u.eval(p) - functionevaluation);
  //    maxErr = (diff > maxErr) ? diff : maxErr;
  //    L2Err += diff * diff;
  //  }
  //  L2Err = sqrt(L2Err / static_cast<double>(numMCpoints));
  //  std::cout << "max err: " << maxErr << " , L2 err: " << L2Err << std::endl;
  //  p = sgpp::base::DataVector(2, 0.4);
  //  std::cout << "u(0.4,0.4) = " << u.eval(p) << std::endl;
  //  std::cout << "|1 - u(0.4,0.4)| = " << fabs(1 - u.eval(p)) << std::endl;
  //----------------------------------------------------------------

  double numGridPoints = static_cast<double>(gridStorage.getSize());

  // calculate mean value via quadrature
  auto quadOperation = createexpUniformBsplineQuadratureCoefficientOperation(degree, numDimensions,
                                                                             coefficientStorage);
  quadOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  double mean = quadOperation->getResult();

  // calculate variance via massMatrix on the SG
  sgpp::base::Grid* gridptr = grid.get();
  sgpp::combigrid::LTwoScalarProductHashMapNakBsplineBoundaryCombigrid massMatrix(gridptr);
  sgpp::base::DataVector product(alpha.size(), 0);
  massMatrix.mult(alpha, product);
  double meanSquare = product.dotProduct(alpha);

  double variance = meanSquare - mean * mean;

  return std::vector<double>{mean, variance, numGridPoints};
}

std::vector<double> evaluateBsplineInterpolant(
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure,
    size_t numDimensions, size_t degree, sgpp::base::DataMatrix params,
    std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage) {
  // create CT interpolation operation
  auto interpolationOperation =
      createBsplineLinearCoefficientOperation(degree, numDimensions, coefficientStorage);

  // convert level structure to SG
  std::shared_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(numDimensions, degree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  convertexpUniformBoundaryCombigridToHierarchicalSparseGrid(levelStructure, gridStorage);

  // interpolate on SG
  sgpp::base::DataMatrix interpolParams(numDimensions, gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    interpolParams.setColumn(i, p);
  }

  // obtain function values from combigrid surrogate
  interpolationOperation->setParameters(interpolParams);
  interpolationOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector f_values = interpolationOperation->getResult();

  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::UMFPACK sleSolver;
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }

  sgpp::optimization::InterpolantScalarFunction u(*grid, alpha);
  sgpp::base::DataVector p(numDimensions);
  std::vector<double> evaluations(params.getNcols(), 0.0);
// ToDo (rehmemk) set NUM_OMP_THREADS according to python. Measure if this omp directive is useful
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < params.getNcols(); i++) {
    params.getColumn(i, p);
    evaluations[i] = u.eval(p);
  }
  return evaluations;
}

sgpp::base::DataVector createInterpolantOnConvertedExpUnifromBoundaryCombigird(
    std::shared_ptr<sgpp::base::Grid>& grid, sgpp::base::GridStorage& gridStorage,
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation>& combigridInterpolationOperation,
    std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> const& levelStructure) {
  sgpp::base::DataMatrix interpolParams(gridStorage.getDimension(), gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    interpolParams.setColumn(i, p);
  }

  // obtain function values from combigrid surrogate
  combigridInterpolationOperation->setParameters(interpolParams);
  //
  combigridInterpolationOperation->getLevelManager()->addLevelsFromStructure(levelStructure);
  sgpp::base::DataVector f_values = combigridInterpolationOperation->getResult();

  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::UMFPACK sleSolver;
  sgpp::base::DataVector alpha(grid->getSize());
  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }
  return alpha;
}

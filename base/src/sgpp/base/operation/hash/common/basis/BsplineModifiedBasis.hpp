// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINE_MODIFIED_BASE_HPP
#define BSPLINE_MODIFIED_BASE_HPP

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

/**
 * Modified B-spline basis on Noboundary grids.
 */
template <class LT, class IT>
class BsplineModifiedBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  BsplineModifiedBasis() : bsplineBasis(BsplineBasis<LT, IT>()) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit BsplineModifiedBasis(size_t degree) : bsplineBasis(BsplineBasis<LT, IT>(degree)) {}

  /**
   * Destructor.
   */
  ~BsplineModifiedBasis() override {}

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @return      value of modified uniform B-spline (e.g. index == 1)
   */
  inline double modifiedBSpline(double x, size_t p) const {
    // wrote the following by hand... nah, not really
    // (thanks to Sage & Python!)
    switch (p) {
      case 3:
        if (x >= 3.0) {
          return 0.0;
        } else if (x < 1.0) {
          return -x + 2.0;
        } else if (x < 2.0) {
          return 1.0 / 6.0 * x * x * x - 0.5 * x * x - 0.5 * x + 11.0 / 6.0;
        } else {
          return -1.0 / 6.0 * x * x * x + 1.5 * x * x - 4.5 * x + 4.5;
        }

        break;

      case 5:

        // use Horner's method for evaluation
        if (x >= 4.0) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 120.0;
          result *= x;
          result *= x;
          result *= x;
          result = -1.0 + result * x;
          result = 2.0 + result * x;
          return result;
        } else if (x < 2.0) {
          double result = -0.025;
          result = 1.0 / 6.0 + result * x;
          result = -1.0 / 3.0 + result * x;
          result = 1.0 / 3.0 + result * x;
          result = -7.0 / 6.0 + result * x;
          result = 61.0 / 30.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 0.025;
          result = -1.0 / 3.0 + result * x;
          result = 5.0 / 3.0 + result * x;
          result = -11.0 / 3.0 + result * x;
          result = 17.0 / 6.0 + result * x;
          result = 13.0 / 30.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 120.0;
          result = 1.0 / 6.0 + result * x;
          result = -4.0 / 3.0 + result * x;
          result = 16.0 / 3.0 + result * x;
          result = -32.0 / 3.0 + result * x;
          result = 128.0 / 15.0 + result * x;
          return result;
        }

        break;

      case 7:
        if (x >= 5.0) {
          return 0.0;
        } else if (x < 1.0) {
          double result = -1.0 / 1008.0;
          result = 1.0 / 720.0 + result * x;
          result = 1.0 / 240.0 + result * x;
          result = 1.0 / 144.0 + result * x;
          result = 1.0 / 144.0 + result * x;
          result = 1.0 / 240.0 + result * x;
          result = -719.0 / 720.0 + result * x;
          result = 10081.0 / 5040.0 + result * x;
          return result;
        } else if (x < 2.0) {
          double result = 1.0 / 504.0;
          result = -7.0 / 360.0 + result * x;
          result = 1.0 / 15.0 + result * x;
          result = -7.0 / 72.0 + result * x;
          result = 1.0 / 9.0 + result * x;
          result = -7.0 / 120.0 + result * x;
          result = -44.0 / 45.0 + result * x;
          result = 719.0 / 360.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = -1.0 / 504.0;
          result = 13.0 / 360.0 + result * x;
          result = -4.0 / 15.0 + result * x;
          result = 73.0 / 72.0 + result * x;
          result = -19.0 / 9.0 + result * x;
          result = 313.0 / 120.0 + result * x;
          result = -124.0 / 45.0 + result * x;
          result = 6313.0 / 2520.0 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = 1.0 / 1008.0;
          result = -19.0 / 720.0 + result * x;
          result = 71.0 / 240.0 + result * x;
          result = -259.0 / 144.0 + result * x;
          result = 911.0 / 144.0 + result * x;
          result = -3019.0 / 240.0 + result * x;
          result = 8951.0 / 720.0 + result * x;
          result = -20179.0 / 5040.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 5040.0;
          result = 1.0 / 144.0 + result * x;
          result = -5.0 / 48.0 + result * x;
          result = 125.0 / 144.0 + result * x;
          result = -625.0 / 144.0 + result * x;
          result = 625.0 / 48.0 + result * x;
          result = -3125.0 / 144.0 + result * x;
          result = 15625.0 / 1008.0 + result * x;
          return result;
        }

        break;

      // for the sake of completeness...
      case 1:
        if (x >= 2.0) {
          return 0.0;
        } else {
          return -x + 2.0;
        }

        break;

      /*case 2:
          if (x >= 2.5)
          {
              return 0.0;
          } else if (x < 1.5)
          {
              return -x + 2.0;
          } else
          {
              return 0.5*x*x - 2.5*x + 3.125;
          }
          break;
      case 4:
          if (x >= 3.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              return -x + 2.0;
          } else if (x < 1.5)
          {
              double result = 1.0/24.0;
              result = -1.0/12.0 + result * x;
              result = 0.0625 + result * x;
              result = -49.0/48.0 + result * x;
              result = 769.0/384.0 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = -1.0/12.0;
              result = 2.0/3.0 + result * x;
              result = -1.625 + result * x;
              result = 2.0/3.0 + result * x;
              result = 263.0/192.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/24.0;
              result = -7.0/12.0 + result * x;
              result = 3.0625 + result * x;
              result = -343.0/48.0 + result * x;
              result = 2401.0/384.0 + result * x;
              return result;
          }
          break;
      case 6:
          if (x >= 4.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              double result = 1.0/720.0;
              result = 1.0/240.0 + result * x;
              result = 1.0/192.0 + result * x;
              result = 1.0/288.0 + result * x;
              result = 1.0/768.0 + result * x;
              result = -3839.0/3840.0 + result * x;
              result = 92161.0/46080.0 + result * x;
              return result;
          } else if (x < 1.5)
          {
              double result = -1.0/180.0;
              result = 0.025 + result * x;
              result = -1.0/48.0 + result * x;
              result = 1.0/48.0 + result * x;
              result = -1.0/192.0 + result * x;
              result = -0.9984375 + result * x;
              result = 23039.0/11520.0 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = 1.0/120.0;
              result = -0.1 + result * x;
              result = 43.0/96.0 + result * x;
              result = -11.0/12.0 + result * x;
              result = 403.0/384.0 + result * x;
              result = -1.63125 + result * x;
              result = 49723.0/23040.0 + result * x;
              return result;
          } else if (x < 3.5)
          {
              double result = -1.0/180.0;
              result = 13.0/120.0 + result * x;
              result = -41.0/48.0 + result * x;
              result = 493.0/144.0 + result * x;
              result = -1361.0/192.0 + result * x;
              result = 12493.0/1920.0 + result * x;
              result = -14201.0/11520.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/720.0;
              result = -0.0375 + result * x;
              result = 0.421875 + result * x;
              result = -2.53125 + result * x;
              result = 8.54296875 + result * x;
              result = -15.37734375 + result * x;
              result = 11.5330078125 + result * x;
              return result;
          }
          break;
      case 8:
          if (x >= 5.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              double result = -1.0/6720.0;
              result = -1.0/2520.0 + result * x;
              result = 1.0/2880.0 + result * x;
              result = 1.0/288.0 + result * x;
              result = 37.0/4608.0 + result * x;
              result = 59.0/5760.0 + result * x;
              result = 361.0/46080.0 + result * x;
              result = -32147.0/32256.0 + result * x;
              result = 10325197.0/5160960.0 + result * x;
              return result;
          } else if (x < 1.5)
          {
              double result = 1.0/2688.0;
              result = -5.0/2016.0 + result * x;
              result = 23.0/5760.0 + result * x;
              result = -1.0/5760.0 + result * x;
              result = 95.0/9216.0 + result * x;
              result = 43.0/4608.0 + result * x;
              result = 743.0/92160.0 + result * x;
              result = -642961.0/645120.0 + result * x;
              result = 4130083.0/2064384.0 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = -1.0/2016.0;
              result = 1.0/126.0 + result * x;
              result = -73.0/1440.0 + result * x;
              result = 59.0/360.0 + result * x;
              result = -685.0/2304.0 + result * x;
              result = 109.0/288.0 + result * x;
              result = -6193.0/23040.0 + result * x;
              result = -35401.0/40320.0 + result * x;
              result = 1021039.0/516096.0 + result * x;
              return result;
          } else if (x < 3.5)
          {
              double result = 1.0/2688.0;
              result = -19.0/2016.0 + result * x;
              result = 583.0/5760.0 + result * x;
              result = -3431.0/5760.0 + result * x;
              result = 19135.0/9216.0 + result * x;
              result = -20131.0/4608.0 + result * x;
              result = 522103.0/92160.0 + result * x;
              result = -3300791.0/645120.0 + result * x;
              result = 6818531.0/2064384.0 + result * x;
              return result;
          } else if (x < 4.5)
          {
              double result = -1.0/6720.0;
              result = 13.0/2520.0 + result * x;
              result = -223.0/2880.0 + result * x;
              result = 943.0/1440.0 + result * x;
              result = -15643.0/4608.0 + result * x;
              result = 63073.0/5760.0 + result * x;
              result = -974263.0/46080.0 + result * x;
              result = 3498403.0/161280.0 + result * x;
              result = -43484083.0/5160960.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/40320.0;
              result = -11.0/10080.0 + result * x;
              result = 121.0/5760.0 + result * x;
              result = -1331.0/5760.0 + result * x;
              result = 14641.0/9216.0 + result * x;
              result = -161051.0/23040.0 + result * x;
              result = 1771561.0/92160.0 + result * x;
              result = -19487171.0/645120.0 + result * x;
              result = 214358881.0/10321920.0 + result * x;
              return result;
          }
          break;
      case 9:
          if (x >= 6.0)
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/17280.0;
              result = -1.0/6720.0 + result * x;
              result = -1.0/2520.0 + result * x;
              result *= x;
              result = 1.0/360.0 + result * x;
              result = 1.0/120.0 + result * x;
              result = 7.0/540.0 + result * x;
              result = 1.0/84.0 + result * x;
              result = -5009.0/5040.0 + result * x;
              result = 1441.0/720.0 + result * x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/10368.0;
              result = 5.0/4032.0 + result * x;
              result = -1.0/168.0 + result * x;
              result = 7.0/540.0 + result * x;
              result = -1.0/60.0 + result * x;
              result = 1.0/36.0 + result * x;
              result *= x;
              result = 11.0/630.0 + result * x;
              result = -209.0/210.0 + result * x;
              result = 1297.0/648.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/10368.0;
              result = -1.0/448.0 + result * x;
              result = 11.0/504.0 + result * x;
              result = -7.0/60.0 + result * x;
              result = 67.0/180.0 + result * x;
              result = -0.75 + result * x;
              result = 28.0/27.0 + result * x;
              result = -61.0/70.0 + result * x;
              result = -347.0/630.0 + result * x;
              result = 137.0/72.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/17280.0;
              result = 13.0/6720.0 + result * x;
              result = -71.0/2520.0 + result * x;
              result = 7.0/30.0 + result * x;
              result = -433.0/360.0 + result * x;
              result = 3.975 + result * x;
              result = -4543.0/540.0 + result * x;
              result = 1579.0/140.0 + result * x;
              result = -48703.0/5040.0 + result * x;
              result = 3557.0/720.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/51840.0;
              result = -17.0/20160.0 + result * x;
              result = 41.0/2520.0 + result * x;
              result = -49.0/270.0 + result * x;
              result = 463.0/360.0 + result * x;
              result = -2153.0/360.0 + result * x;
              result = 9793.0/540.0 + result * x;
              result = -43133.0/1260.0 + result * x;
              result = 180673.0/5040.0 + result * x;
              result = -99059.0/6480.0 + result * x;
              return result;
          } else
          {
              double result = -1.0/362880.0;
              result = 1.0/6720.0 + result * x;
              result = -1.0/280.0 + result * x;
              result = 0.05 + result * x;
              result = -0.45 + result * x;
              result = 2.7 + result * x;
              result = -10.8 + result * x;
              result = 972.0/35.0 + result * x;
              result = -1458.0/35.0 + result * x;
              result = 972.0/35.0 + result * x;
              return result;
          }
          break;*/
      default:
        // the degree is too damn high
        // ==> calculate modified B-spline value via definition
        double y = 0.0;
        double x2 = x + static_cast<double>(p + 1) / 2.0 - 1.0;

        if (x2 > static_cast<double>(p) + 1.0) {
          return 0.0;
        }

        // the upper summation bound is defined to be ceil((p + 1) / 2.0),
        // which is the same as (p + 2) / 2 written in C
        for (size_t k = 0; k <= (p + 2) / 2; k++) {
          // x2 is chosen such that it holds: x2 = x + (p+1)/2 + k - 1
          y += static_cast<double>(k + 1) * bsplineBasis.uniformBSpline(x2, p);
          // the rounding errors induced by this method can be neglected
          x2 += 1.0;
        }

        return y;
    }
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @return      value of derivative of modified uniform B-spline
   *              (e.g. index == 1)
   */
  inline double modifiedBSplineDx(double x, size_t p) const {
    switch (p) {
      case 3:
        if (x >= 3.0) {
          return 0.0;
        } else if (x < 1.0) {
          return -1.0;
        } else if (x < 2.0) {
          return 0.5 * x * x - x - 0.5;
        } else {
          return -0.5 * x * x + 3.0 * x - 4.5;
        }

        break;

      case 5:
        if (x >= 4.0) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 24.0;
          result *= x;
          result *= x;
          result *= x;
          result = -1.0 + result * x;
          return result;
        } else if (x < 2.0) {
          double result = -0.125;
          result = 2.0 / 3.0 + result * x;
          result = -1.0 + result * x;
          result = 2.0 / 3.0 + result * x;
          result = -7.0 / 6.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 0.125;
          result = -4.0 / 3.0 + result * x;
          result = 5.0 + result * x;
          result = -22.0 / 3.0 + result * x;
          result = 17.0 / 6.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 24.0;
          result = 2.0 / 3.0 + result * x;
          result = -4.0 + result * x;
          result = 32.0 / 3.0 + result * x;
          result = -32.0 / 3.0 + result * x;
          return result;
        }

        break;

      case 7:
        if (x >= 5.0) {
          return 0.0;
        } else if (x < 1.0) {
          double result = -1.0 / 144.0;
          result = 1.0 / 120.0 + result * x;
          result = 1.0 / 48.0 + result * x;
          result = 1.0 / 36.0 + result * x;
          result = 1.0 / 48.0 + result * x;
          result = 1.0 / 120.0 + result * x;
          result = -719.0 / 720.0 + result * x;
          return result;
        } else if (x < 2.0) {
          double result = 1.0 / 72.0;
          result = -7.0 / 60.0 + result * x;
          result = 1.0 / 3.0 + result * x;
          result = -7.0 / 18.0 + result * x;
          result = 1.0 / 3.0 + result * x;
          result = -7.0 / 60.0 + result * x;
          result = -44.0 / 45.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = -1.0 / 72.0;
          result = 13.0 / 60.0 + result * x;
          result = -4.0 / 3.0 + result * x;
          result = 73.0 / 18.0 + result * x;
          result = -19.0 / 3.0 + result * x;
          result = 313.0 / 60.0 + result * x;
          result = -124.0 / 45.0 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = 1.0 / 144.0;
          result = -19.0 / 120.0 + result * x;
          result = 71.0 / 48.0 + result * x;
          result = -259.0 / 36.0 + result * x;
          result = 911.0 / 48.0 + result * x;
          result = -3019.0 / 120.0 + result * x;
          result = 8951.0 / 720.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 720.0;
          result = 1.0 / 24.0 + result * x;
          result = -25.0 / 48.0 + result * x;
          result = 125.0 / 36.0 + result * x;
          result = -625.0 / 48.0 + result * x;
          result = 625.0 / 24.0 + result * x;
          result = -3125.0 / 144.0 + result * x;
          return result;
        }

        break;

      case 1:
        if (x >= 2.0) {
          return 0.0;
        } else {
          return -1.0;
        }

        break;

      /*case 2:
          if (x >= 2.5)
          {
              return 0.0;
          } else if (x < 1.5)
          {
              return -1.0;
          } else
          {
              return x - 2.5;
          }
          break;
      case 4:
          if (x >= 3.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              return -1.0;
          } else if (x < 1.5)
          {
              return 1.0/6.0*x*x*x - 0.25*x*x + 0.125*x - 49.0/48.0;
          } else if (x < 2.5)
          {
              return -1.0/3.0*x*x*x + 2.0*x*x - 3.25*x + 2.0/3.0;
          } else
          {
              return 1.0/6.0*x*x*x - 1.75*x*x + 6.125*x - 343.0/48.0;
          }
          break;
      case 6:
          if (x >= 4.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              double result = 1.0/120.0;
              result = 1.0/48.0 + result * x;
              result = 1.0/48.0 + result * x;
              result = 1.0/96.0 + result * x;
              result = 1.0/384.0 + result * x;
              result = -3839.0/3840.0 + result * x;
              return result;
          } else if (x < 1.5)
          {
              double result = -1.0/30.0;
              result = 0.125 + result * x;
              result = -1.0/12.0 + result * x;
              result = 0.0625 + result * x;
              result = -1.0/96.0 + result * x;
              result = -0.9984375 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = 0.05;
              result = -0.5 + result * x;
              result = 43.0/24.0 + result * x;
              result = -2.75 + result * x;
              result = 403.0/192.0 + result * x;
              result = -1.63125 + result * x;
              return result;
          } else if (x < 3.5)
          {
              double result = -1.0/30.0;
              result = 13.0/24.0 + result * x;
              result = -41.0/12.0 + result * x;
              result = 493.0/48.0 + result * x;
              result = -1361.0/96.0 + result * x;
              result = 12493.0/1920.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/120.0;
              result = -0.1875 + result * x;
              result = 1.6875 + result * x;
              result = -7.59375 + result * x;
              result = 17.0859375 + result * x;
              result = -15.37734375 + result * x;
              return result;
          }
          break;
      case 8:
          if (x >= 5.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              double result = -1.0/840.0;
              result = -1.0/360.0 + result * x;
              result = 1.0/480.0 + result * x;
              result = 5.0/288.0 + result * x;
              result = 37.0/1152.0 + result * x;
              result = 59.0/1920.0 + result * x;
              result = 361.0/23040.0 + result * x;
              result = -32147.0/32256.0 + result * x;
              return result;
          } else if (x < 1.5)
          {
              double result = 1.0/336.0;
              result = -5.0/288.0 + result * x;
              result = 23.0/960.0 + result * x;
              result = -1.0/1152.0 + result * x;
              result = 95.0/2304.0 + result * x;
              result = 43.0/1536.0 + result * x;
              result = 743.0/46080.0 + result * x;
              result = -642961.0/645120.0 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = -1.0/252.0;
              result = 1.0/18.0 + result * x;
              result = -73.0/240.0 + result * x;
              result = 59.0/72.0 + result * x;
              result = -685.0/576.0 + result * x;
              result = 109.0/96.0 + result * x;
              result = -6193.0/11520.0 + result * x;
              result = -35401.0/40320.0 + result * x;
              return result;
          } else if (x < 3.5)
          {
              double result = 1.0/336.0;
              result = -19.0/288.0 + result * x;
              result = 583.0/960.0 + result * x;
              result = -3431.0/1152.0 + result * x;
              result = 19135.0/2304.0 + result * x;
              result = -20131.0/1536.0 + result * x;
              result = 522103.0/46080.0 + result * x;
              result = -3300791.0/645120.0 + result * x;
              return result;
          } else if (x < 4.5)
          {
              double result = -1.0/840.0;
              result = 13.0/360.0 + result * x;
              result = -223.0/480.0 + result * x;
              result = 943.0/288.0 + result * x;
              result = -15643.0/1152.0 + result * x;
              result = 63073.0/1920.0 + result * x;
              result = -974263.0/23040.0 + result * x;
              result = 3498403.0/161280.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/5040.0;
              result = -11.0/1440.0 + result * x;
              result = 121.0/960.0 + result * x;
              result = -1331.0/1152.0 + result * x;
              result = 14641.0/2304.0 + result * x;
              result = -161051.0/7680.0 + result * x;
              result = 1771561.0/46080.0 + result * x;
              result = -19487171.0/645120.0 + result * x;
              return result;
          }
          break;
      case 9:
          if (x >= 6.0)
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/1920.0;
              result = -1.0/840.0 + result * x;
              result = -1.0/360.0 + result * x;
              result *= x;
              result = 1.0/72.0 + result * x;
              result = 1.0/30.0 + result * x;
              result = 7.0/180.0 + result * x;
              result = 1.0/42.0 + result * x;
              result = -5009.0/5040.0 + result * x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/1152.0;
              result = 5.0/504.0 + result * x;
              result = -1.0/24.0 + result * x;
              result = 7.0/90.0 + result * x;
              result = -1.0/12.0 + result * x;
              result = 1.0/9.0 + result * x;
              result *= x;
              result = 11.0/315.0 + result * x;
              result = -209.0/210.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/1152.0;
              result = -1.0/56.0 + result * x;
              result = 11.0/72.0 + result * x;
              result = -0.7 + result * x;
              result = 67.0/36.0 + result * x;
              result = -3.0 + result * x;
              result = 28.0/9.0 + result * x;
              result = -61.0/35.0 + result * x;
              result = -347.0/630.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/1920.0;
              result = 13.0/840.0 + result * x;
              result = -71.0/360.0 + result * x;
              result = 1.4 + result * x;
              result = -433.0/72.0 + result * x;
              result = 15.9 + result * x;
              result = -4543.0/180.0 + result * x;
              result = 1579.0/70.0 + result * x;
              result = -48703.0/5040.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/5760.0;
              result = -17.0/2520.0 + result * x;
              result = 41.0/360.0 + result * x;
              result = -49.0/45.0 + result * x;
              result = 463.0/72.0 + result * x;
              result = -2153.0/90.0 + result * x;
              result = 9793.0/180.0 + result * x;
              result = -43133.0/630.0 + result * x;
              result = 180673.0/5040.0 + result * x;
              return result;
          } else
          {
              double result = -1.0/40320.0;
              result = 1.0/840.0 + result * x;
              result = -0.025 + result * x;
              result = 0.3 + result * x;
              result = -2.25 + result * x;
              result = 10.8 + result * x;
              result = -32.4 + result * x;
              result = 1944.0/35.0 + result * x;
              result = -1458.0/35.0 + result * x;
              return result;
          }
          break;*/
      default:
        double y = 0.0;
        double x2 = x + static_cast<double>(p + 1) / 2.0 - 1.0;

        if (x2 > static_cast<double>(p) + 1.0) {
          return 0.0;
        }

        for (size_t k = 0; k <= (p + 2) / 2; k++) {
          y += static_cast<double>(k + 1) * bsplineBasis.uniformBSplineDx(x2, p);
          x2 += 1.0;
        }

        return y;
    }
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @return      value of 2nd derivative of modified uniform B-spline
   *              (e.g. index == 1)
   */
  inline double modifiedBSplineDxDx(double x, size_t p) const {
    switch (p) {
      case 3:
        if (x >= 3.0) {
          return 0.0;
        } else if (x < 1.0) {
          return 0.0;
        } else if (x < 2.0) {
          return x - 1.0;
        } else {
          return -x + 3.0;
        }

        break;

      case 5:
        if (x >= 4.0) {
          return 0.0;
        } else if (x < 1.0) {
          return 1.0 / 6.0 * x * x * x;
        } else if (x < 2.0) {
          return -0.5 * x * x * x + 2.0 * x * x - 2.0 * x + 2.0 / 3.0;
        } else if (x < 3.0) {
          return 0.5 * x * x * x - 4.0 * x * x + 10.0 * x - 22.0 / 3.0;
        } else {
          return -1.0 / 6.0 * x * x * x + 2.0 * x * x - 8.0 * x + 32.0 / 3.0;
        }

        break;

      case 7:
        if (x >= 5.0) {
          return 0.0;
        } else if (x < 1.0) {
          double result = -1.0 / 24.0;
          result = 1.0 / 24.0 + result * x;
          result = 1.0 / 12.0 + result * x;
          result = 1.0 / 12.0 + result * x;
          result = 1.0 / 24.0 + result * x;
          result = 1.0 / 120.0 + result * x;
          return result;
        } else if (x < 2.0) {
          double result = 1.0 / 12.0;
          result = -7.0 / 12.0 + result * x;
          result = 4.0 / 3.0 + result * x;
          result = -7.0 / 6.0 + result * x;
          result = 2.0 / 3.0 + result * x;
          result = -7.0 / 60.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = -1.0 / 12.0;
          result = 13.0 / 12.0 + result * x;
          result = -16.0 / 3.0 + result * x;
          result = 73.0 / 6.0 + result * x;
          result = -38.0 / 3.0 + result * x;
          result = 313.0 / 60.0 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = 1.0 / 24.0;
          result = -19.0 / 24.0 + result * x;
          result = 71.0 / 12.0 + result * x;
          result = -259.0 / 12.0 + result * x;
          result = 911.0 / 24.0 + result * x;
          result = -3019.0 / 120.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 120.0;
          result = 5.0 / 24.0 + result * x;
          result = -25.0 / 12.0 + result * x;
          result = 125.0 / 12.0 + result * x;
          result = -625.0 / 24.0 + result * x;
          result = 625.0 / 24.0 + result * x;
          return result;
        }

        break;

      case 1:
        return 0.0;
        break;

      /*case 2:
          if (x >= 2.5)
          {
              return 0.0;
          } else if (x < 1.5)
          {
              return 0.0;
          } else
          {
              return 1.0;
          }
          break;
      case 4:
          if ((x < 0.0) || (x >= 3.5))
          {
              return 0.0;
          } else if (x < 0.5)
          {
              return 0.0;
          } else if (x < 1.5)
          {
              return 0.5*x*x - 0.5*x + 0.125;
          } else if (x < 2.5)
          {
              return -x*x + 4.0*x - 3.25;
          } else
          {
              return 0.5*x*x - 3.5*x + 6.125;
          }
          break;
      case 6:
          if (x >= 4.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              double result = 1.0/24.0;
              result = 1.0/12.0 + result * x;
              result = 0.0625 + result * x;
              result = 1.0/48.0 + result * x;
              result = 1.0/384.0 + result * x;
              return result;
          } else if (x < 1.5)
          {
              double result = -1.0/6.0;
              result = 0.5 + result * x;
              result = -0.25 + result * x;
              result = 0.125 + result * x;
              result = -1.0/96.0 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = 0.25;
              result = -2.0 + result * x;
              result = 5.375 + result * x;
              result = -5.5 + result * x;
              result = 403.0/192.0 + result * x;
              return result;
          } else if (x < 3.5)
          {
              double result = -1.0/6.0;
              result = 13.0/6.0 + result * x;
              result = -10.25 + result * x;
              result = 493.0/24.0 + result * x;
              result = -1361.0/96.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/24.0;
              result = -0.75 + result * x;
              result = 5.0625 + result * x;
              result = -15.1875 + result * x;
              result = 17.0859375 + result * x;
              return result;
          }
          break;
      case 8:
          if (x >= 5.5)
          {
              return 0.0;
          } else if (x < 0.5)
          {
              double result = -1.0/120.0;
              result = -1.0/60.0 + result * x;
              result = 1.0/96.0 + result * x;
              result = 5.0/72.0 + result * x;
              result = 37.0/384.0 + result * x;
              result = 59.0/960.0 + result * x;
              result = 361.0/23040.0 + result * x;
              return result;
          } else if (x < 1.5)
          {
              double result = 1.0/48.0;
              result = -5.0/48.0 + result * x;
              result = 23.0/192.0 + result * x;
              result = -1.0/288.0 + result * x;
              result = 95.0/768.0 + result * x;
              result = 43.0/768.0 + result * x;
              result = 743.0/46080.0 + result * x;
              return result;
          } else if (x < 2.5)
          {
              double result = -1.0/36.0;
              result = 1.0/3.0 + result * x;
              result = -73.0/48.0 + result * x;
              result = 59.0/18.0 + result * x;
              result = -685.0/192.0 + result * x;
              result = 109.0/48.0 + result * x;
              result = -6193.0/11520.0 + result * x;
              return result;
          } else if (x < 3.5)
          {
              double result = 1.0/48.0;
              result = -19.0/48.0 + result * x;
              result = 583.0/192.0 + result * x;
              result = -3431.0/288.0 + result * x;
              result = 19135.0/768.0 + result * x;
              result = -20131.0/768.0 + result * x;
              result = 522103.0/46080.0 + result * x;
              return result;
          } else if (x < 4.5)
          {
              double result = -1.0/120.0;
              result = 13.0/60.0 + result * x;
              result = -223.0/96.0 + result * x;
              result = 943.0/72.0 + result * x;
              result = -15643.0/384.0 + result * x;
              result = 63073.0/960.0 + result * x;
              result = -974263.0/23040.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/720.0;
              result = -11.0/240.0 + result * x;
              result = 121.0/192.0 + result * x;
              result = -1331.0/288.0 + result * x;
              result = 14641.0/768.0 + result * x;
              result = -161051.0/3840.0 + result * x;
              result = 1771561.0/46080.0 + result * x;
              return result;
          }
          break;
      case 9:
          if (x >= 6.0)
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/240.0;
              result = -1.0/120.0 + result * x;
              result = -1.0/60.0 + result * x;
              result *= x;
              result = 1.0/18.0 + result * x;
              result = 0.1 + result * x;
              result = 7.0/90.0 + result * x;
              result = 1.0/42.0 + result * x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/144.0;
              result = 5.0/72.0 + result * x;
              result = -0.25 + result * x;
              result = 7.0/18.0 + result * x;
              result = -1.0/3.0 + result * x;
              result = 1.0/3.0 + result * x;
              result *= x;
              result = 11.0/315.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/144.0;
              result = -0.125 + result * x;
              result = 11.0/12.0 + result * x;
              result = -3.5 + result * x;
              result = 67.0/9.0 + result * x;
              result = -9.0 + result * x;
              result = 56.0/9.0 + result * x;
              result = -61.0/35.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/240.0;
              result = 13.0/120.0 + result * x;
              result = -71.0/60.0 + result * x;
              result = 7.0 + result * x;
              result = -433.0/18.0 + result * x;
              result = 47.7 + result * x;
              result = -4543.0/90.0 + result * x;
              result = 1579.0/70.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/720.0;
              result = -17.0/360.0 + result * x;
              result = 41.0/60.0 + result * x;
              result = -49.0/9.0 + result * x;
              result = 463.0/18.0 + result * x;
              result = -2153.0/30.0 + result * x;
              result = 9793.0/90.0 + result * x;
              result = -43133.0/630.0 + result * x;
              return result;
          } else
          {
              double result = -1.0/5040.0;
              result = 1.0/120.0 + result * x;
              result = -0.15 + result * x;
              result = 1.5 + result * x;
              result = -9.0 + result * x;
              result = 32.4 + result * x;
              result = -64.8 + result * x;
              result = 1944.0/35.0 + result * x;
              return result;
          }
          break;*/
      default:
        double y = 0.0;
        double x2 = x + static_cast<double>(p + 1) / 2.0 - 1.0;

        if (x2 > static_cast<double>(p) + 1.0) {
          return 0.0;
        }

        for (size_t k = 0; k <= (p + 2) / 2; k++) {
          y += static_cast<double>(k + 1) * bsplineBasis.uniformBSplineDxDx(x2, p);
          x2 += 1.0;
        }

        return y;
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of modified B-spline basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 1) {
      return 1.0;
    }

    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);

    if (i == hInv - 1) {
      // mirror the situation at x = 0.5
      x = 1.0 - x;
      i = 1;
    }

    if (i == 1) {
      return modifiedBSpline(x * hInvDbl, bsplineBasis.getDegree());
    } else {
      return bsplineBasis.uniformBSpline(
          x * hInvDbl - static_cast<double>(i) +
              static_cast<double>(bsplineBasis.getDegree() + 1) / 2.0,
          bsplineBasis.getDegree());
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of modified
   *              B-spline basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    if (l == 1) {
      return 0.0;
    }

    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    // inner derivative
    double dxFactor = hInvDbl;

    if (i == hInv - 1) {
      // mirror the situation at x = 0.5
      // (don't forget the inner derivative!)
      x = 1.0 - x;
      i = 1;
      dxFactor *= -1.0;
    }

    if (i == 1) {
      return dxFactor * modifiedBSplineDx(x * hInvDbl, bsplineBasis.getDegree());
    } else {
      return dxFactor *
             bsplineBasis.uniformBSplineDx(
                 x * hInvDbl - static_cast<double>(i) +
                     static_cast<double>(bsplineBasis.getDegree() + 1) / 2.0,
                 bsplineBasis.getDegree());
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of modified
   *              B-spline basis function
   */
  inline double evalDxDx(LT l, IT i, double x) {
    if (l == 1) {
      return 0.0;
    }

    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    // inner derivative
    double dxFactor = hInvDbl;

    if (i == hInv - 1) {
      // mirror the situation at x = 0.5
      x = 1.0 - x;
      i = 1;
    }

    if (i == 1) {
      return dxFactor * dxFactor * modifiedBSplineDxDx(x * hInvDbl, bsplineBasis.getDegree());
    } else {
      return dxFactor * dxFactor *
             bsplineBasis.uniformBSplineDxDx(
                 x * hInvDbl - static_cast<double>(i) +
                     static_cast<double>(bsplineBasis.getDegree() + 1) / 2.0,
                 bsplineBasis.getDegree());
    }
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return bsplineBasis.getDegree(); }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @return      value of the integral
   */
  inline double getIntegral(LT l, IT i) override {
    if (l == 1) {
      return 1;
    }
    const IT hInv = static_cast<IT>(1) << l;
    const double hInvDbl = static_cast<double>(hInv);
    if (i == 1 || i == hInv - 1) {
      switch (bsplineBasis.getDegree()) {
        case 1:
          return 2 / hInvDbl;
        case 3:
          return (25.0 / 12.0) / hInvDbl;
        case 5:
          // 1081/720 + 0.581944 + 0.0819444 + 1/720
          return (2.0 + 1.0 / 6.0) / hInvDbl;
        case 7:
          // 20243/13440 + 2495/4032 + 479/4032 + 83/13440 + 1/40320 (last part is cut off on level
          // 2)
          if (l == 2) {
            return 90718.0 / 40320.0 / hInvDbl;  // = 2.25 - 2/40320
          } else {
            return 90719.0 / 40320.0 / hInvDbl;  // = 2.25 - 1/4032
          }
        default:
          throw operation_exception(
              "BsplineModifiedBasis::getIntegral() "
              "only implemented for 1 <= degree <= 7");
          break;
      }
    } else {
      return bsplineBasis.getIntegral(l, i);
    }
  }

 protected:
  /// B-spline basis for B-spline evaluation
  BsplineBasis<LT, IT> bsplineBasis;
};

// default type-def (unsigned int for level and index)
typedef BsplineModifiedBasis<unsigned int, unsigned int> SBsplineModifiedBase;

}  // namespace base
}  // namespace sgpp

#endif /* BSPLINE_MODIFIED_BASE_HPP */

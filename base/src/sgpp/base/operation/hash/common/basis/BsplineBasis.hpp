// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BSPLINE_BASE_HPP
#define BSPLINE_BASE_HPP

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <iostream>

namespace sgpp {
namespace base {

/**
 * B-spline basis on Noboundary grids.
 */
template <class LT, class IT>
class BsplineBasis : public Basis<LT, IT> {
 public:
  /**
   * Default constructor.
   */
  BsplineBasis() : degree(0) {}

  /**
   * Constructor.
   *
   * @param degree    B-spline degree, must be odd
   *                  (if it's even, degree - 1 is used)
   */
  explicit BsplineBasis(size_t degree) : degree(degree) {
    if (degree < 1) {
      this->degree = 1;
    } else if (degree % 2 == 0) {
      this->degree = degree - 1;
    }
  }

  /**
   * Destructor.
   */
  ~BsplineBasis() override {}

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @return      value of uniform B-spline
   *              (with knots \f$\{0, 1, ..., p+1\}\f$)
   */
  inline double uniformBSpline(double x, size_t p) const {
    // wrote the following by hand... nah, not really
    // (thanks to Sage & Python!)
    switch (p) {
      case 3:
        if ((x < 0.0) || (x >= 4.0)) {
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

      case 5:

        // use Horner's method for evaluation
        if ((x < 0.0) || (x >= 6.0)) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 120.0;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          return result;
        } else if (x < 2.0) {
          double result = -1.0 / 24.0;
          result = 0.25 + result * x;
          result = -0.5 + result * x;
          result = 0.5 + result * x;
          result = -0.25 + result * x;
          result = 0.05 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 1.0 / 12.0;
          result = -1.0 + result * x;
          result = 4.5 + result * x;
          result = -9.5 + result * x;
          result = 9.75 + result * x;
          result = -3.95 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = -1.0 / 12.0;
          result = 1.5 + result * x;
          result = -10.5 + result * x;
          result = 35.5 + result * x;
          result = -57.75 + result * x;
          result = 36.55 + result * x;
          return result;
        } else if (x < 5.0) {
          double result = 1.0 / 24.0;
          result = -1.0 + result * x;
          result = 9.5 + result * x;
          result = -44.5 + result * x;
          result = 102.25 + result * x;
          result = -91.45 + result * x;
          return result;
        } else {
          double result = -1.0 / 120.0;
          result = 0.25 + result * x;
          result = -3.0 + result * x;
          result = 18.0 + result * x;
          result = -54.0 + result * x;
          result = 64.8 + result * x;
          return result;
        }

        break;

      case 7:
        if ((x < 0.0) || (x >= 8.0)) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 5040.0;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          return result;
        } else if (x < 2.0) {
          double result = -1.0 / 720.0;
          result = 1.0 / 90.0 + result * x;
          result = -1.0 / 30.0 + result * x;
          result = 1.0 / 18.0 + result * x;
          result = -1.0 / 18.0 + result * x;
          result = 1.0 / 30.0 + result * x;
          result = -1.0 / 90.0 + result * x;
          result = 1.0 / 630.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 1.0 / 240.0;
          result = -1.0 / 15.0 + result * x;
          result = 13.0 / 30.0 + result * x;
          result = -1.5 + result * x;
          result = 55.0 / 18.0 + result * x;
          result = -3.7 + result * x;
          result = 223.0 / 90.0 + result * x;
          result = -149.0 / 210.0 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = -1.0 / 144.0;
          result = 1.0 / 6.0 + result * x;
          result = -5.0 / 3.0 + result * x;
          result = 9.0 + result * x;
          result = -256.0 / 9.0 + result * x;
          result = 53.0 + result * x;
          result = -488.0 / 9.0 + result * x;
          result = 2477.0 / 105.0 + result * x;
          return result;
        } else if (x < 5.0) {
          double result = 1.0 / 144.0;
          result = -2.0 / 9.0 + result * x;
          result = 3.0 + result * x;
          result = -199.0 / 9.0 + result * x;
          result = 96.0 + result * x;
          result = -737.0 / 3.0 + result * x;
          result = 344.0 + result * x;
          result = -64249.0 / 315.0 + result * x;
          return result;
        } else if (x < 6.0) {
          double result = -1.0 / 240.0;
          result = 1.0 / 6.0 + result * x;
          result = -17.0 / 6.0 + result * x;
          result = 26.5 + result * x;
          result = -2647.0 / 18.0 + result * x;
          result = 483.5 + result * x;
          result = -15683.0 / 18.0 + result * x;
          result = 139459.0 / 210.0 + result * x;
          return result;
        } else if (x < 7.0) {
          double result = 1.0 / 720.0;
          result = -1.0 / 15.0 + result * x;
          result = 41.0 / 30.0 + result * x;
          result = -15.5 + result * x;
          result = 1889.0 / 18.0 + result * x;
          result = -423.7 + result * x;
          result = 84881.0 / 90.0 + result * x;
          result = -187133.0 / 210.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 5040.0;
          result = 1.0 / 90.0 + result * x;
          result = -4.0 / 15.0 + result * x;
          result = 32.0 / 9.0 + result * x;
          result = -256.0 / 9.0 + result * x;
          result = 2048.0 / 15.0 + result * x;
          result = -16384.0 / 45.0 + result * x;
          result = 131072.0 / 315.0 + result * x;
          return result;
        }

        break;

      case 0:
        if ((x < 0.0) || (x >= 1.0)) {
          return 0.0;
        } else {
          return 1.0;
        }

        break;

      case 1:
        if ((x < 0.0) || (x >= 2.0)) {
          return 0.0;
        } else if (x < 1.0) {
          return x;
        } else {
          return -x + 2.0;
        }

        break;

      // for the sake of completeness...
      /*case 2:
          if ((x < 0.0) || (x >= 3.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              return 0.5*x*x;
          } else if (x < 2.0)
          {
              return -x*x + 3.0*x - 1.5;
          } else
          {
              return 0.5*x*x - 3.0*x + 4.5;
          }
          break;
      case 4:
          if ((x < 0.0) || (x >= 5.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/24.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/6.0;
              result = 5.0/6.0 + result * x;
              result = -1.25 + result * x;
              result = 5.0/6.0 + result * x;
              result = -5.0/24.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 0.25;
              result = -2.5 + result * x;
              result = 8.75 + result * x;
              result = -12.5 + result * x;
              result = 155.0/24.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/6.0;
              result = 2.5 + result * x;
              result = -13.75 + result * x;
              result = 32.5 + result * x;
              result = -655.0/24.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/24.0;
              result = -5.0/6.0 + result * x;
              result = 6.25 + result * x;
              result = -125.0/6.0 + result * x;
              result = 625.0/24.0 + result * x;
              return result;
          }
          break;
      case 6:
          if ((x < 0.0) || (x >= 7.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/720.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/120.0;
              result = 7.0/120.0 + result * x;
              result = -7.0/48.0 + result * x;
              result = 7.0/36.0 + result * x;
              result = -7.0/48.0 + result * x;
              result = 7.0/120.0 + result * x;
              result = -7.0/720.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/48.0;
              result = -7.0/24.0 + result * x;
              result = 77.0/48.0 + result * x;
              result = -161.0/36.0 + result * x;
              result = 329.0/48.0 + result * x;
              result = -133.0/24.0 + result * x;
              result = 1337.0/720.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/36.0;
              result = 7.0/12.0 + result * x;
              result = -119.0/24.0 + result * x;
              result = 196.0/9.0 + result * x;
              result = -1253.0/24.0 + result * x;
              result = 196.0/3.0 + result * x;
              result = -12089.0/360.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/48.0;
              result = -7.0/12.0 + result * x;
              result = 161.0/24.0 + result * x;
              result = -364.0/9.0 + result * x;
              result = 3227.0/24.0 + result * x;
              result = -700.0/3.0 + result * x;
              result = 59591.0/360.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -1.0/120.0;
              result = 7.0/24.0 + result * x;
              result = -203.0/48.0 + result * x;
              result = 1169.0/36.0 + result * x;
              result = -6671.0/48.0 + result * x;
              result = 7525.0/24.0 + result * x;
              result = -208943.0/720.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/720.0;
              result = -7.0/120.0 + result * x;
              result = 49.0/48.0 + result * x;
              result = -343.0/36.0 + result * x;
              result = 2401.0/48.0 + result * x;
              result = -16807.0/120.0 + result * x;
              result = 117649.0/720.0 + result * x;
              return result;
          }
          break;
      case 8:
          if ((x < 0.0) || (x >= 9.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/40320.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/5040.0;
              result = 1.0/560.0 + result * x;
              result = -0.00625 + result * x;
              result = 0.0125 + result * x;
              result = -0.015625 + result * x;
              result = 0.0125 + result * x;
              result = -0.00625 + result * x;
              result = 1.0/560.0 + result * x;
              result = -1.0/4480.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/1440.0;
              result = -0.0125 + result * x;
              result = 0.09375 + result * x;
              result = -0.3875 + result * x;
              result = 0.984375 + result * x;
              result = -1.5875 + result * x;
              result = 1.59375 + result * x;
              result = -0.9125 + result * x;
              result = 1023.0/4480.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/720.0;
              result = 0.0375 + result * x;
              result = -0.43125 + result * x;
              result = 2.7625 + result * x;
              result = -10.828125 + result * x;
              result = 26.7625 + result * x;
              result = -40.93125 + result * x;
              result = 35.5375 + result * x;
              result = -60213.0/4480.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/576.0;
              result = -0.0625 + result * x;
              result = 0.96875 + result * x;
              result = -8.4375 + result * x;
              result = 45.171875 + result * x;
              result = -152.4375 + result * x;
              result = 317.46875 + result * x;
              result = -374.0625 + result * x;
              result = 857291.0/4480.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -1.0/720.0;
              result = 0.0625 + result * x;
              result = -1.21875 + result * x;
              result = 13.4375 + result * x;
              result = -91.546875 + result * x;
              result = 394.4375 + result * x;
              result = -1049.71875 + result * x;
              result = 1579.0625 + result * x;
              result = -4611459.0/4480.0 + result * x;
              return result;
          } else if (x < 7.0)
          {
              double result = 1.0/1440.0;
              result = -0.0375 + result * x;
              result = 0.88125 + result * x;
              result = -11.7625 + result * x;
              result = 97.453125 + result * x;
              result = -512.7625 + result * x;
              result = 1671.88125 + result * x;
              result = -3086.5375 + result * x;
              result = 11064957.0/4480.0 + result * x;
              return result;
          } else if (x < 8.0)
          {
              double result = -1.0/5040.0;
              result = 0.0125 + result * x;
              result = -0.34375 + result * x;
              result = 5.3875 + result * x;
              result = -52.609375 + result * x;
              result = 327.5875 + result * x;
              result = -1269.34375 + result * x;
              result = 2795.9125 + result * x;
              result = -11994247.0/4480.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/40320.0;
              result = -1.0/560.0 + result * x;
              result = 0.05625 + result * x;
              result = -1.0125 + result * x;
              result = 11.390625 + result * x;
              result = -82.0125 + result * x;
              result = 369.05625 + result * x;
              result = -531441.0/560.0 + result * x;
              result = 4782969.0/4480.0 + result * x;
              return result;
          }
          break;
      case 9:
          if ((x < 0.0) || (x >= 10.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/362880.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/40320.0;
              result = 1.0/4032.0 + result * x;
              result = -1.0/1008.0 + result * x;
              result = 1.0/432.0 + result * x;
              result = -1.0/288.0 + result * x;
              result = 1.0/288.0 + result * x;
              result = -1.0/432.0 + result * x;
              result = 1.0/1008.0 + result * x;
              result = -1.0/4032.0 + result * x;
              result = 1.0/36288.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/10080.0;
              result = -1.0/504.0 + result * x;
              result = 17.0/1008.0 + result * x;
              result = -35.0/432.0 + result * x;
              result = 71.0/288.0 + result * x;
              result = -143.0/288.0 + result * x;
              result = 287.0/432.0 + result * x;
              result = -575.0/1008.0 + result * x;
              result = 1151.0/4032.0 + result * x;
              result = -329.0/5184.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/4320.0;
              result = 1.0/144.0 + result * x;
              result = -13.0/144.0 + result * x;
              result = 289.0/432.0 + result * x;
              result = -901.0/288.0 + result * x;
              result = 2773.0/288.0 + result * x;
              result = -8461.0/432.0 + result * x;
              result = 3667.0/144.0 + result * x;
              result = -11083.0/576.0 + result * x;
              result = 233893.0/36288.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/2880.0;
              result = -1.0/72.0 + result * x;
              result = 35.0/144.0 + result * x;
              result = -1055.0/432.0 + result * x;
              result = 4475.0/288.0 + result * x;
              result = -18731.0/288.0 + result * x;
              result = 77555.0/432.0 + result * x;
              result = -45485.0/144.0 + result * x;
              result = 185525.0/576.0 + result * x;
              result = -5271131.0/36288.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -1.0/2880.0;
              result = 5.0/288.0 + result * x;
              result = -55.0/144.0 + result * x;
              result = 2095.0/432.0 + result * x;
              result = -11275.0/288.0 + result * x;
              result = 60019.0/288.0 + result * x;
              result = -316195.0/432.0 + result * x;
              result = 235765.0/144.0 + result * x;
              result = -1220725.0/576.0 + result * x;
              result = 43947619.0/36288.0 + result * x;
              return result;
          } else if (x < 7.0)
          {
              double result = 1.0/4320.0;
              result = -1.0/72.0 + result * x;
              result = 53.0/144.0 + result * x;
              result = -2441.0/432.0 + result * x;
              result = 15941.0/288.0 + result * x;
              result = -103277.0/288.0 + result * x;
              result = 663581.0/432.0 + result * x;
              result = -604043.0/144.0 + result * x;
              result = 3818123.0/576.0 + result * x;
              result = -167683997.0/36288.0 + result * x;
              return result;
          } else if (x < 8.0)
          {
              double result = -1.0/10080.0;
              result = 1.0/144.0 + result * x;
              result = -31.0/144.0 + result * x;
              result = 1675.0/432.0 + result * x;
              result = -12871.0/288.0 + result * x;
              result = 98407.0/288.0 + result * x;
              result = -748207.0/432.0 + result * x;
              result = 807745.0/144.0 + result * x;
              result = -6064393.0/576.0 + result * x;
              result = 316559287.0/36288.0 + result * x;
              return result;
          } else if (x < 9.0)
          {
              double result = 1.0/40320.0;
              result = -1.0/504.0 + result * x;
              result = 71.0/1008.0 + result * x;
              result = -629.0/432.0 + result * x;
              result = 5561.0/288.0 + result * x;
              result = -49049.0/288.0 + result * x;
              result = 431441.0/432.0 + result * x;
              result = -3782969.0/1008.0 + result * x;
              result = 33046721.0/4032.0 + result * x;
              result = -287420489.0/36288.0 + result * x;
              return result;
          } else
          {
              double result = -1.0/362880.0;
              result = 1.0/4032.0 + result * x;
              result = -5.0/504.0 + result * x;
              result = 25.0/108.0 + result * x;
              result = -125.0/36.0 + result * x;
              result = 625.0/18.0 + result * x;
              result = -6250.0/27.0 + result * x;
              result = 62500.0/63.0 + result * x;
              result = -156250.0/63.0 + result * x;
              result = 1562500.0/567.0 + result * x;
              return result;
          }
          break;*/
      default:
        // the degree is too damn high
        // ==> calculate B-spline value by Cox-de-Boor recursion
        const double pDbl = static_cast<double>(p);

        if ((x < 0.0) || (x >= pDbl + 1.0)) {
          return 0.0;
        } else {
          return (x / pDbl) * uniformBSpline(x, p - 1) +
                 ((pDbl + 1.0 - x) / pDbl) * uniformBSpline(x - 1.0, p - 1);
        }
    }
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @return      value of derivative of uniform B-spline
   *              (with knots \f$\{0, 1, ..., p+1\}\f$)
   */
  inline double uniformBSplineDx(double x, size_t p) const {
    switch (p) {
      case 3:
        if ((x < 0.0) || (x >= 4.0)) {
          return 0.0;
        } else if (x < 1.0) {
          return 0.5 * x * x;
        } else if (x < 2.0) {
          return -1.5 * x * x + 4.0 * x - 2.0;
        } else if (x < 3.0) {
          return 1.5 * x * x - 8.0 * x + 10.0;
        } else {
          return -0.5 * x * x + 4.0 * x - 8.0;
        }

        break;

      case 5:
        if ((x < 0.0) || (x >= 6.0)) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 24.0;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          return result;
        } else if (x < 2.0) {
          double result = -5.0 / 24.0;
          result = 1.0 + result * x;
          result = -1.5 + result * x;
          result = 1.0 + result * x;
          result = -0.25 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 5.0 / 12.0;
          result = -4.0 + result * x;
          result = 13.5 + result * x;
          result = -19.0 + result * x;
          result = 9.75 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = -5.0 / 12.0;
          result = 6.0 + result * x;
          result = -31.5 + result * x;
          result = 71.0 + result * x;
          result = -57.75 + result * x;
          return result;
        } else if (x < 5.0) {
          double result = 5.0 / 24.0;
          result = -4.0 + result * x;
          result = 28.5 + result * x;
          result = -89.0 + result * x;
          result = 102.25 + result * x;
          return result;
        } else {
          double result = -1.0 / 24.0;
          result = 1.0 + result * x;
          result = -9.0 + result * x;
          result = 36.0 + result * x;
          result = -54.0 + result * x;
          return result;
        }

        break;

      case 7:
        if ((x < 0.0) || (x >= 8.0)) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 720.0;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          return result;
        } else if (x < 2.0) {
          double result = -7.0 / 720.0;
          result = 1.0 / 15.0 + result * x;
          result = -1.0 / 6.0 + result * x;
          result = 2.0 / 9.0 + result * x;
          result = -1.0 / 6.0 + result * x;
          result = 1.0 / 15.0 + result * x;
          result = -1.0 / 90.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 7.0 / 240.0;
          result = -0.4 + result * x;
          result = 13.0 / 6.0 + result * x;
          result = -6.0 + result * x;
          result = 55.0 / 6.0 + result * x;
          result = -7.4 + result * x;
          result = 223.0 / 90.0 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = -7.0 / 144.0;
          result = 1.0 + result * x;
          result = -25.0 / 3.0 + result * x;
          result = 36.0 + result * x;
          result = -256.0 / 3.0 + result * x;
          result = 106.0 + result * x;
          result = -488.0 / 9.0 + result * x;
          return result;
        } else if (x < 5.0) {
          double result = 7.0 / 144.0;
          result = -4.0 / 3.0 + result * x;
          result = 15.0 + result * x;
          result = -796.0 / 9.0 + result * x;
          result = 288.0 + result * x;
          result = -1474.0 / 3.0 + result * x;
          result = 344.0 + result * x;
          return result;
        } else if (x < 6.0) {
          double result = -7.0 / 240.0;
          result = 1.0 + result * x;
          result = -85.0 / 6.0 + result * x;
          result = 106.0 + result * x;
          result = -2647.0 / 6.0 + result * x;
          result = 967.0 + result * x;
          result = -15683.0 / 18.0 + result * x;
          return result;
        } else if (x < 7.0) {
          double result = 7.0 / 720.0;
          result = -0.4 + result * x;
          result = 41.0 / 6.0 + result * x;
          result = -62.0 + result * x;
          result = 1889.0 / 6.0 + result * x;
          result = -847.4 + result * x;
          result = 84881.0 / 90.0 + result * x;
          return result;
        } else {
          double result = -1.0 / 720.0;
          result = 1.0 / 15.0 + result * x;
          result = -4.0 / 3.0 + result * x;
          result = 128.0 / 9.0 + result * x;
          result = -256.0 / 3.0 + result * x;
          result = 4096.0 / 15.0 + result * x;
          result = -16384.0 / 45.0 + result * x;
          return result;
        }

        break;

      case 0:
        return 0.0;
        break;

      case 1:
        if ((x < 0.0) || (x >= 2.0)) {
          return 0.0;
        } else if (x < 1.0) {
          return 1.0;
        } else {
          return -1.0;
        }

        break;

      /*case 2:
          if ((x < 0.0) || (x >= 3.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              return x;
          } else if (x < 2.0)
          {
              return -2.0*x + 3.0;
          } else
          {
              return x - 3.0;
          }
          break;
      case 4:
          if ((x < 0.0) || (x >= 5.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              return 1.0/6.0*x*x*x;
          } else if (x < 2.0)
          {
              return -2.0/3.0*x*x*x + 2.5*x*x - 2.5*x + 5.0/6.0;
          } else if (x < 3.0)
          {
              return x*x*x - 7.5*x*x + 17.5*x - 12.5;
          } else if (x < 4.0)
          {
              return -2.0/3.0*x*x*x + 7.5*x*x - 27.5*x + 32.5;
          } else
          {
              return 1.0/6.0*x*x*x - 2.5*x*x + 12.5*x - 125.0/6.0;
          }
          break;
      case 6:
          if ((x < 0.0) || (x >= 7.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/120.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -0.05;
              result = 7.0/24.0 + result * x;
              result = -7.0/12.0 + result * x;
              result = 7.0/12.0 + result * x;
              result = -7.0/24.0 + result * x;
              result = 7.0/120.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 0.125;
              result = -35.0/24.0 + result * x;
              result = 77.0/12.0 + result * x;
              result = -161.0/12.0 + result * x;
              result = 329.0/24.0 + result * x;
              result = -133.0/24.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/6.0;
              result = 35.0/12.0 + result * x;
              result = -119.0/6.0 + result * x;
              result = 196.0/3.0 + result * x;
              result = -1253.0/12.0 + result * x;
              result = 196.0/3.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 0.125;
              result = -35.0/12.0 + result * x;
              result = 161.0/6.0 + result * x;
              result = -364.0/3.0 + result * x;
              result = 3227.0/12.0 + result * x;
              result = -700.0/3.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -0.05;
              result = 35.0/24.0 + result * x;
              result = -203.0/12.0 + result * x;
              result = 1169.0/12.0 + result * x;
              result = -6671.0/24.0 + result * x;
              result = 7525.0/24.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/120.0;
              result = -7.0/24.0 + result * x;
              result = 49.0/12.0 + result * x;
              result = -343.0/12.0 + result * x;
              result = 2401.0/24.0 + result * x;
              result = -16807.0/120.0 + result * x;
              return result;
          }
          break;
      case 8:
          if ((x < 0.0) || (x >= 9.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/5040.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/630.0;
              result = 0.0125 + result * x;
              result = -0.0375 + result * x;
              result = 0.0625 + result * x;
              result = -0.0625 + result * x;
              result = 0.0375 + result * x;
              result = -0.0125 + result * x;
              result = 1.0/560.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/180.0;
              result = -0.0875 + result * x;
              result = 0.5625 + result * x;
              result = -1.9375 + result * x;
              result = 3.9375 + result * x;
              result = -4.7625 + result * x;
              result = 3.1875 + result * x;
              result = -0.9125 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/90.0;
              result = 0.2625 + result * x;
              result = -2.5875 + result * x;
              result = 13.8125 + result * x;
              result = -43.3125 + result * x;
              result = 80.2875 + result * x;
              result = -81.8625 + result * x;
              result = 35.5375 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 1.0/72.0;
              result = -0.4375 + result * x;
              result = 5.8125 + result * x;
              result = -42.1875 + result * x;
              result = 180.6875 + result * x;
              result = -457.3125 + result * x;
              result = 634.9375 + result * x;
              result = -374.0625 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -1.0/90.0;
              result = 0.4375 + result * x;
              result = -7.3125 + result * x;
              result = 67.1875 + result * x;
              result = -366.1875 + result * x;
              result = 1183.3125 + result * x;
              result = -2099.4375 + result * x;
              result = 1579.0625 + result * x;
              return result;
          } else if (x < 7.0)
          {
              double result = 1.0/180.0;
              result = -0.2625 + result * x;
              result = 5.2875 + result * x;
              result = -58.8125 + result * x;
              result = 389.8125 + result * x;
              result = -1538.2875 + result * x;
              result = 3343.7625 + result * x;
              result = -3086.5375 + result * x;
              return result;
          } else if (x < 8.0)
          {
              double result = -1.0/630.0;
              result = 0.0875 + result * x;
              result = -2.0625 + result * x;
              result = 26.9375 + result * x;
              result = -210.4375 + result * x;
              result = 982.7625 + result * x;
              result = -2538.6875 + result * x;
              result = 2795.9125 + result * x;
              return result;
          } else
          {
              double result = 1.0/5040.0;
              result = -0.0125 + result * x;
              result = 0.3375 + result * x;
              result = -5.0625 + result * x;
              result = 45.5625 + result * x;
              result = -246.0375 + result * x;
              result = 738.1125 + result * x;
              result = -531441.0/560.0 + result * x;
              return result;
          }
          break;
      case 9:
          if ((x < 0.0) || (x >= 10.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/40320.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/4480.0;
              result = 1.0/504.0 + result * x;
              result = -1.0/144.0 + result * x;
              result = 1.0/72.0 + result * x;
              result = -5.0/288.0 + result * x;
              result = 1.0/72.0 + result * x;
              result = -1.0/144.0 + result * x;
              result = 1.0/504.0 + result * x;
              result = -1.0/4032.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/1120.0;
              result = -1.0/63.0 + result * x;
              result = 17.0/144.0 + result * x;
              result = -35.0/72.0 + result * x;
              result = 355.0/288.0 + result * x;
              result = -143.0/72.0 + result * x;
              result = 287.0/144.0 + result * x;
              result = -575.0/504.0 + result * x;
              result = 1151.0/4032.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/480.0;
              result = 1.0/18.0 + result * x;
              result = -91.0/144.0 + result * x;
              result = 289.0/72.0 + result * x;
              result = -4505.0/288.0 + result * x;
              result = 2773.0/72.0 + result * x;
              result = -8461.0/144.0 + result * x;
              result = 3667.0/72.0 + result * x;
              result = -11083.0/576.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 0.003125;
              result = -1.0/9.0 + result * x;
              result = 245.0/144.0 + result * x;
              result = -1055.0/72.0 + result * x;
              result = 22375.0/288.0 + result * x;
              result = -18731.0/72.0 + result * x;
              result = 77555.0/144.0 + result * x;
              result = -45485.0/72.0 + result * x;
              result = 185525.0/576.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -0.003125;
              result = 5.0/36.0 + result * x;
              result = -385.0/144.0 + result * x;
              result = 2095.0/72.0 + result * x;
              result = -56375.0/288.0 + result * x;
              result = 60019.0/72.0 + result * x;
              result = -316195.0/144.0 + result * x;
              result = 235765.0/72.0 + result * x;
              result = -1220725.0/576.0 + result * x;
              return result;
          } else if (x < 7.0)
          {
              double result = 1.0/480.0;
              result = -1.0/9.0 + result * x;
              result = 371.0/144.0 + result * x;
              result = -2441.0/72.0 + result * x;
              result = 79705.0/288.0 + result * x;
              result = -103277.0/72.0 + result * x;
              result = 663581.0/144.0 + result * x;
              result = -604043.0/72.0 + result * x;
              result = 3818123.0/576.0 + result * x;
              return result;
          } else if (x < 8.0)
          {
              double result = -1.0/1120.0;
              result = 1.0/18.0 + result * x;
              result = -217.0/144.0 + result * x;
              result = 1675.0/72.0 + result * x;
              result = -64355.0/288.0 + result * x;
              result = 98407.0/72.0 + result * x;
              result = -748207.0/144.0 + result * x;
              result = 807745.0/72.0 + result * x;
              result = -6064393.0/576.0 + result * x;
              return result;
          } else if (x < 9.0)
          {
              double result = 1.0/4480.0;
              result = -1.0/63.0 + result * x;
              result = 71.0/144.0 + result * x;
              result = -629.0/72.0 + result * x;
              result = 27805.0/288.0 + result * x;
              result = -49049.0/72.0 + result * x;
              result = 431441.0/144.0 + result * x;
              result = -3782969.0/504.0 + result * x;
              result = 33046721.0/4032.0 + result * x;
              return result;
          } else
          {
              double result = -1.0/40320.0;
              result = 1.0/504.0 + result * x;
              result = -5.0/72.0 + result * x;
              result = 25.0/18.0 + result * x;
              result = -625.0/36.0 + result * x;
              result = 1250.0/9.0 + result * x;
              result = -6250.0/9.0 + result * x;
              result = 125000.0/63.0 + result * x;
              result = -156250.0/63.0 + result * x;
              return result;
          }
          break;*/
      default:
        if ((x < 0.0) || (x >= static_cast<double>(p) + 1.0)) {
          return 0.0;
        } else {
          return uniformBSpline(x, p - 1) - uniformBSpline(x - 1.0, p - 1);
        }
    }
  }

  /**
   * @param x     evaluation point
   * @param p     B-spline degree
   * @return      value of 2nd derivative of uniform B-spline
   *              (with knots \f$\{0, 1, ..., p+1\}\f$)
   */
  inline double uniformBSplineDxDx(double x, size_t p) const {
    switch (p) {
      case 3:
        if ((x < 0.0) || (x >= 4.0)) {
          return 0.0;
        } else if (x < 1.0) {
          return x;
        } else if (x < 2.0) {
          return -3.0 * x + 4.0;
        } else if (x < 3.0) {
          return 3.0 * x - 8.0;
        } else {
          return -x + 4.0;
        }

        break;

      case 5:
        if ((x < 0.0) || (x >= 6.0)) {
          return 0.0;
        } else if (x < 1.0) {
          return 1.0 / 6.0 * x * x * x;
        } else if (x < 2.0) {
          return -5.0 / 6.0 * x * x * x + 3.0 * x * x - 3.0 * x + 1.0;
        } else if (x < 3.0) {
          return 5.0 / 3.0 * x * x * x - 12.0 * x * x + 27.0 * x - 19.0;
        } else if (x < 4.0) {
          return -5.0 / 3.0 * x * x * x + 18.0 * x * x - 63.0 * x + 71.0;
        } else if (x < 5.0) {
          return 5.0 / 6.0 * x * x * x - 12.0 * x * x + 57.0 * x - 89.0;
        } else {
          return -1.0 / 6.0 * x * x * x + 3.0 * x * x - 18.0 * x + 36.0;
        }

        break;

      case 7:
        if ((x < 0.0) || (x >= 8.0)) {
          return 0.0;
        } else if (x < 1.0) {
          double result = 1.0 / 120.0;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          result *= x;
          return result;
        } else if (x < 2.0) {
          double result = -7.0 / 120.0;
          result = 1.0 / 3.0 + result * x;
          result = -2.0 / 3.0 + result * x;
          result = 2.0 / 3.0 + result * x;
          result = -1.0 / 3.0 + result * x;
          result = 1.0 / 15.0 + result * x;
          return result;
        } else if (x < 3.0) {
          double result = 0.175;
          result = -2.0 + result * x;
          result = 26.0 / 3.0 + result * x;
          result = -18.0 + result * x;
          result = 55.0 / 3.0 + result * x;
          result = -7.4 + result * x;
          return result;
        } else if (x < 4.0) {
          double result = -7.0 / 24.0;
          result = 5.0 + result * x;
          result = -100.0 / 3.0 + result * x;
          result = 108.0 + result * x;
          result = -512.0 / 3.0 + result * x;
          result = 106.0 + result * x;
          return result;
        } else if (x < 5.0) {
          double result = 7.0 / 24.0;
          result = -20.0 / 3.0 + result * x;
          result = 60.0 + result * x;
          result = -796.0 / 3.0 + result * x;
          result = 576.0 + result * x;
          result = -1474.0 / 3.0 + result * x;
          return result;
        } else if (x < 6.0) {
          double result = -0.175;
          result = 5.0 + result * x;
          result = -170.0 / 3.0 + result * x;
          result = 318.0 + result * x;
          result = -2647.0 / 3.0 + result * x;
          result = 967.0 + result * x;
          return result;
        } else if (x < 7.0) {
          double result = 7.0 / 120.0;
          result = -2.0 + result * x;
          result = 82.0 / 3.0 + result * x;
          result = -186.0 + result * x;
          result = 1889.0 / 3.0 + result * x;
          result = -847.4 + result * x;
          return result;
        } else {
          double result = -1.0 / 120.0;
          result = 1.0 / 3.0 + result * x;
          result = -16.0 / 3.0 + result * x;
          result = 128.0 / 3.0 + result * x;
          result = -512.0 / 3.0 + result * x;
          result = 4096.0 / 15.0 + result * x;
          return result;
        }

        break;

      case 0:
      case 1:
        return 0.0;
        break;

      /*case 2:
          if ((x < 0.0) || (x >= 3.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              return 1.0;
          } else if (x < 2.0)
          {
              return -2.0;
          } else
          {
              return 1.0;
          }
          break;
      case 4:
          if ((x < 0.0) || (x >= 5.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              return 0.5*x*x;
          } else if (x < 2.0)
          {
              return -2.0*x*x + 5.0*x - 2.5;
          } else if (x < 3.0)
          {
              return 3.0*x*x - 15.0*x + 17.5;
          } else if (x < 4.0)
          {
              return -2.0*x*x + 15.0*x - 27.5;
          } else
          {
              return 0.5*x*x - 5.0*x + 12.5;
          }
          break;
      case 6:
          if ((x < 0.0) || (x >= 7.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/24.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -0.25;
              result = 7.0/6.0 + result * x;
              result = -1.75 + result * x;
              result = 7.0/6.0 + result * x;
              result = -7.0/24.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 0.625;
              result = -35.0/6.0 + result * x;
              result = 19.25 + result * x;
              result = -161.0/6.0 + result * x;
              result = 329.0/24.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -5.0/6.0;
              result = 35.0/3.0 + result * x;
              result = -59.5 + result * x;
              result = 392.0/3.0 + result * x;
              result = -1253.0/12.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 0.625;
              result = -35.0/3.0 + result * x;
              result = 80.5 + result * x;
              result = -728.0/3.0 + result * x;
              result = 3227.0/12.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -0.25;
              result = 35.0/6.0 + result * x;
              result = -50.75 + result * x;
              result = 1169.0/6.0 + result * x;
              result = -6671.0/24.0 + result * x;
              return result;
          } else
          {
              double result = 1.0/24.0;
              result = -7.0/6.0 + result * x;
              result = 12.25 + result * x;
              result = -343.0/6.0 + result * x;
              result = 2401.0/24.0 + result * x;
              return result;
          }
          break;
      case 8:
          if ((x < 0.0) || (x >= 9.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/720.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/90.0;
              result = 0.075 + result * x;
              result = -0.1875 + result * x;
              result = 0.25 + result * x;
              result = -0.1875 + result * x;
              result = 0.075 + result * x;
              result = -0.0125 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 7.0/180.0;
              result = -0.525 + result * x;
              result = 2.8125 + result * x;
              result = -7.75 + result * x;
              result = 11.8125 + result * x;
              result = -9.525 + result * x;
              result = 3.1875 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -7.0/90.0;
              result = 1.575 + result * x;
              result = -12.9375 + result * x;
              result = 55.25 + result * x;
              result = -129.9375 + result * x;
              result = 160.575 + result * x;
              result = -81.8625 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 7.0/72.0;
              result = -2.625 + result * x;
              result = 29.0625 + result * x;
              result = -168.75 + result * x;
              result = 542.0625 + result * x;
              result = -914.625 + result * x;
              result = 634.9375 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -7.0/90.0;
              result = 2.625 + result * x;
              result = -36.5625 + result * x;
              result = 268.75 + result * x;
              result = -1098.5625 + result * x;
              result = 2366.625 + result * x;
              result = -2099.4375 + result * x;
              return result;
          } else if (x < 7.0)
          {
              double result = 7.0/180.0;
              result = -1.575 + result * x;
              result = 26.4375 + result * x;
              result = -235.25 + result * x;
              result = 1169.4375 + result * x;
              result = -3076.575 + result * x;
              result = 3343.7625 + result * x;
              return result;
          } else if (x < 8.0)
          {
              double result = -1.0/90.0;
              result = 0.525 + result * x;
              result = -10.3125 + result * x;
              result = 107.75 + result * x;
              result = -631.3125 + result * x;
              result = 1965.525 + result * x;
              result = -2538.6875 + result * x;
              return result;
          } else
          {
              double result = 1.0/720.0;
              result = -0.075 + result * x;
              result = 1.6875 + result * x;
              result = -20.25 + result * x;
              result = 136.6875 + result * x;
              result = -492.075 + result * x;
              result = 738.1125 + result * x;
              return result;
          }
          break;
      case 9:
          if ((x < 0.0) || (x >= 10.0))
          {
              return 0.0;
          } else if (x < 1.0)
          {
              double result = 1.0/5040.0;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              result *= x;
              return result;
          } else if (x < 2.0)
          {
              double result = -1.0/560.0;
              result = 1.0/72.0 + result * x;
              result = -1.0/24.0 + result * x;
              result = 5.0/72.0 + result * x;
              result = -5.0/72.0 + result * x;
              result = 1.0/24.0 + result * x;
              result = -1.0/72.0 + result * x;
              result = 1.0/504.0 + result * x;
              return result;
          } else if (x < 3.0)
          {
              double result = 1.0/140.0;
              result = -1.0/9.0 + result * x;
              result = 17.0/24.0 + result * x;
              result = -175.0/72.0 + result * x;
              result = 355.0/72.0 + result * x;
              result = -143.0/24.0 + result * x;
              result = 287.0/72.0 + result * x;
              result = -575.0/504.0 + result * x;
              return result;
          } else if (x < 4.0)
          {
              double result = -1.0/60.0;
              result = 7.0/18.0 + result * x;
              result = -91.0/24.0 + result * x;
              result = 1445.0/72.0 + result * x;
              result = -4505.0/72.0 + result * x;
              result = 2773.0/24.0 + result * x;
              result = -8461.0/72.0 + result * x;
              result = 3667.0/72.0 + result * x;
              return result;
          } else if (x < 5.0)
          {
              double result = 0.025;
              result = -7.0/9.0 + result * x;
              result = 245.0/24.0 + result * x;
              result = -5275.0/72.0 + result * x;
              result = 22375.0/72.0 + result * x;
              result = -18731.0/24.0 + result * x;
              result = 77555.0/72.0 + result * x;
              result = -45485.0/72.0 + result * x;
              return result;
          } else if (x < 6.0)
          {
              double result = -0.025;
              result = 35.0/36.0 + result * x;
              result = -385.0/24.0 + result * x;
              result = 10475.0/72.0 + result * x;
              result = -56375.0/72.0 + result * x;
              result = 60019.0/24.0 + result * x;
              result = -316195.0/72.0 + result * x;
              result = 235765.0/72.0 + result * x;
              return result;
          } else if (x < 7.0)
          {
              double result = 1.0/60.0;
              result = -7.0/9.0 + result * x;
              result = 371.0/24.0 + result * x;
              result = -12205.0/72.0 + result * x;
              result = 79705.0/72.0 + result * x;
              result = -103277.0/24.0 + result * x;
              result = 663581.0/72.0 + result * x;
              result = -604043.0/72.0 + result * x;
              return result;
          } else if (x < 8.0)
          {
              double result = -1.0/140.0;
              result = 7.0/18.0 + result * x;
              result = -217.0/24.0 + result * x;
              result = 8375.0/72.0 + result * x;
              result = -64355.0/72.0 + result * x;
              result = 98407.0/24.0 + result * x;
              result = -748207.0/72.0 + result * x;
              result = 807745.0/72.0 + result * x;
              return result;
          } else if (x < 9.0)
          {
              double result = 1.0/560.0;
              result = -1.0/9.0 + result * x;
              result = 71.0/24.0 + result * x;
              result = -3145.0/72.0 + result * x;
              result = 27805.0/72.0 + result * x;
              result = -49049.0/24.0 + result * x;
              result = 431441.0/72.0 + result * x;
              result = -3782969.0/504.0 + result * x;
              return result;
          } else
          {
              double result = -1.0/5040.0;
              result = 1.0/72.0 + result * x;
              result = -5.0/12.0 + result * x;
              result = 125.0/18.0 + result * x;
              result = -625.0/9.0 + result * x;
              result = 1250.0/3.0 + result * x;
              result = -12500.0/9.0 + result * x;
              result = 125000.0/63.0 + result * x;
              return result;
          }
          break;*/
      default:
        if ((x < 0.0) || (x >= static_cast<double>(p) + 1.0)) {
          return 0.0;
        } else {
          return uniformBSpline(x, p - 2) - 2.0 * uniformBSpline(x - 1.0, p - 2) +
                 uniformBSpline(x - 2.0, p - 2);
        }
    }
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of B-spline basis function
   */
  inline double eval(LT l, IT i, double x) override {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);

    return uniformBSpline(
        x * hInv - static_cast<double>(i) + static_cast<double>(this->degree + 1) / 2.0,
        this->degree);
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of derivative of B-spline basis function
   */
  inline double evalDx(LT l, IT i, double x) override {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);

    return hInv * uniformBSplineDx(x * hInv - static_cast<double>(i) +
                                       static_cast<double>(this->degree + 1) / 2.0,
                                   this->degree);
  }

  /**
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of 2nd derivative of B-spline basis function
   */
  inline double evalDxDx(LT l, IT i, double x) {
    const double hInv = static_cast<double>(static_cast<IT>(1) << l);

    return hInv * hInv *
           uniformBSplineDxDx(
               x * hInv - static_cast<double>(i) + static_cast<double>(this->degree + 1) / 2.0,
               this->degree);
  }

  /**
   * @return      B-spline degree
   */
  inline size_t getDegree() const override { return degree; }

  /**
   * @param level     level of basis function
   * @param index     index of basis function
   * @return          value of the Integral
   */
  inline double getIntegral(LT level, IT index) override {
    size_t erster_abschnitt = std::max(0, -static_cast<int>(index - (degree + 1) / 2));
    size_t letzter_abschnitt = std::min(degree, (1 << level) + (degree + 1) / 2 - index - 1);
    double res = 0.0;
    for (size_t j = erster_abschnitt; j <= letzter_abschnitt; j++) {
      switch (degree) {
        case 1:
          switch (j) {
            case 0:
            case 1:
              res += 0.5;
              break;
            default:
              break;
          }
          break;
        case 3:
          switch (j) {
            case 0:
            case 3:
              res += 1.0 / 24.0;
              break;
            case 1:
            case 2:
              res += 11.0 / 24.0;
              break;
            default:
              break;
          }
          break;
        case 5:
          switch (j) {
            case 0:
            case 5:
              res += 1.0 / 720.0;
              break;
            case 1:
            case 4:
              res += 57.0 / 720.0;
              break;
            case 2:
            case 3:
              res += 302.0 / 720.0;
              break;
            default:
              break;
          }
          break;
        case 7:
          switch (j) {
            case 0:
            case 7:
              res += 1.0 / 40320.0;
              break;
            case 1:
            case 6:
              res += 247.0 / 40320.0;
              break;
            case 2:
            case 5:
              res += 4293.0 / 40320.0;
              break;
            case 3:
            case 4:
              res += 15619.0 / 40320.0;
              break;
            default:
              break;
          }
          break;
        default:
          throw operation_exception(
              "BsplineBasis::getIntegral() only"
              "implemented for 1 <= degree <= 7");
          break;
      }
    }
    return 1.0 / (1 << level) * res;
  }

 protected:
  /// degree of the B-spline
  size_t degree;
};

// default type-def (unsigned int for level and index)
typedef BsplineBasis<unsigned int, unsigned int> SBsplineBase;

}  // namespace base
}  // namespace sgpp

#endif /* BSPLINE_BASE_HPP */

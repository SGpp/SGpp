// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>

#include <boost/test/unit_test.hpp>

#include <limits>
#include <utility>
#include <vector>

#include "BasisEval.hpp"

using sgpp::base::BoundingBox1D;
using sgpp::base::DataVector;
using sgpp::base::GridPoint;
using sgpp::base::GridStorage;
using sgpp::base::index_t;
using sgpp::base::level_t;
using sgpp::base::SBasis;
using sgpp::base::Stretching;
using sgpp::base::Stretching1D;

void basisTest(SBasis& basis, const std::vector<level_t>& levels,
               const std::vector<index_t>& indices,
               const std::vector<double>& points,
               const std::vector<double>& testValues) {
  const size_t n = indices.size();
  BOOST_CHECK_EQUAL(levels.size(), n);
  BOOST_CHECK_EQUAL(points.size(), n);
  BOOST_CHECK_EQUAL(testValues.size(), n);

  for (size_t i = 0; i < n; i++) {
    double val = basis.eval(levels[i], indices[i], points[i]);
    BOOST_CHECK_CLOSE(val, testValues[i], 1e-10);
  }
}

void stretchedBasisTest(sgpp::base::SLinearStretchedBase& basis,
                        const std::vector<double>& p,
                        const std::vector<double>& pos0,
                        const std::vector<double>& pos1,
                        const std::vector<double>& testValues) {
  const size_t n = p.size();
  BOOST_CHECK_EQUAL(pos0.size(), n);
  BOOST_CHECK_EQUAL(pos1.size(), n);
  BOOST_CHECK_EQUAL(testValues.size(), n);

  for (size_t i = 0; i < n; i++) {
    double val = basis.stretchedEval(p[i], pos0[i], pos1[i]);
    BOOST_CHECK_CLOSE(val, testValues[i], 1e-10);
  }
}

void linearUniformUnmodifiedTest(SBasis& basis) {
  const std::vector<level_t> levels = {1, 1, 1, 1, 1, 2, 2, 3, 3, 3};
  const std::vector<index_t> indices = {1, 1, 1, 1, 1, 1, 3, 1, 1, 1};

  const std::vector<double> points = {0.5,  0.75, 0.875, 0.0,   1.0,
                                      0.75, 0.75, 0.0,   0.125, 0.25};

  const std::vector<double> testValues = {1.0, 0.5, 0.25, 0.0, 0.0,
                                          0.0, 1.0, 0.0,  1.0, 0.0};

  basisTest(basis, levels, indices, points, testValues);
}

void linearModifiedTest(SBasis& basis) {
  const std::vector<level_t> levels = {1, 1, 1, 1, 1, 2, 2, 2,
                                       2, 2, 2, 2, 3, 3, 3};
  const std::vector<index_t> indices = {1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 3, 3, 1, 1, 1};

  const std::vector<double> points = {0.5,  0.75,  0.875, 0.0,   1.0,
                                      0.0,  0.125, 0.25,  0.375, 0.75,
                                      0.75, 1.0,   0.0,   0.125, 0.25};

  const std::vector<double> testValues = {1.0, 1.0, 1.0, 1.0, 1.0,
                                          2.0, 1.5, 1.0, 0.5, 0.0,
                                          1.0, 2.0, 2.0, 1.0, 0.0};

  basisTest(basis, levels, indices, points, testValues);
}

void linearLevelZeroTest(SBasis& basis) {
  // Test boundary linear basis functions (level 0).
  for (index_t i = 0; i <= 1; i++) {
    for (size_t j = 0; j <= 4; j++) {
      const double x = static_cast<double>(j) / 4.0;
      BOOST_CHECK_EQUAL(
          basis.eval(0, i, x),
          static_cast<double>(i) * x + static_cast<double>(1 - i) * (1.0 - x));
    }
  }
}

double ccKnot(level_t l, index_t i) {
  // Return Clenshaw-Curtis knot with given level and index.
  return 0.5 * (std::cos(M_PI * (1.0 -
                                 static_cast<double>(i) /
                                     std::pow(2.0, static_cast<double>(l)))) +
                1.0);
}

void linearClenshawCurtisBoundaryTest(SBasis& basis) {
  const std::vector<level_t> levels = {1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3};
  const std::vector<level_t> indices = {1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 1, 1};

  const std::vector<double> points = {0.5, 0.75,  0.875,        0.0,
                                      1.0, 0.75,  0.75,         ccKnot(2, 3),
                                      0.0, 0.125, ccKnot(3, 1), 0.25};

  const std::vector<double> testValuesDouble = {
      1.0, 0.5, 0.25, 0.0, 0.0, 0.0,
      0.25 / (static_cast<double>(ccKnot(2, 3)) - 0.5), 1.0, 0.0,
      1.0 -
          (0.125 - static_cast<double>(ccKnot(3, 1))) /
              (static_cast<double>(ccKnot(3, 2) - ccKnot(3, 1))),
      1.0, 0.0};

  const std::vector<double> testValues(testValuesDouble.begin(),
                                       testValuesDouble.end());

  basisTest(basis, levels, indices, points, testValues);
}

void errorTest(double x, double y, double tol) {
  if (std::abs(x) >= 10.0) {
    BOOST_CHECK_SMALL((x - y) / x, tol);
  } else {
    BOOST_CHECK_SMALL(x - y, 10.0 * tol);
  }
}

void derivativesTest(SBasis& basis, size_t degree = 2, level_t startLevel = 1,
                     size_t maxDiscontinuitiesCount = 0) {
  // Test derivatives (up to order deg, max. 2) of basis functions for
  // level >= start_level. Allow for max_discontinuities_count
  // discontinuities (e.g. for Wavelets which are cut off).
  const double dx = 1e-8;
  const double tol1 = 1e-3;
  const double tol2 = 1e-2;
  const double discontinuityTol = 1e5;

  // levels
  for (level_t l = 0; l < 6; l++) {
    // indices
    for (index_t i = 1; i < (static_cast<index_t>(1) << l); i += 2) {
      if (degree >= 1) {
        // test first derivative at boundary (central difference quotient)
        errorTest(
            (basis.eval(l, i, 2.0 * dx) - basis.eval(l, i, 0.0)) / (2.0 * dx),
            basisEvalDx(basis, l, i, dx), tol1);
        errorTest((basis.eval(l, i, 1.0) - basis.eval(l, i, 1.0 - 2.0 * dx)) /
                      (2.0 * dx),
                  basisEvalDx(basis, l, i, 1.0 - dx), tol1);
      }

      if (degree >= 2) {
        // test second derivative at boundary (central difference quotient)
        errorTest((basisEvalDx(basis, l, i, 2.0 * dx) -
                   basisEvalDx(basis, l, i, 0.0)) /
                      (2.0 * dx),
                  basisEvalDxDx(basis, l, i, dx), tol2);
        errorTest((basisEvalDx(basis, l, i, 1.0) -
                   basisEvalDx(basis, l, i, 1.0 - 2.0 * dx)) /
                      (2.0 * dx),
                  basisEvalDxDx(basis, l, i, 1.0 - dx), tol2);
      }

      size_t discontinuities = 0;

      for (size_t j = 1; j < 100; j++) {
        const double x = static_cast<double>(j) / 100.0;

        if (std::abs(basis.eval(l, i, x + dx) - basis.eval(l, i, x - dx)) >
            discontinuityTol * dx) {
          // discontinuity found
          discontinuities++;
        } else {
          // test derivatives only if the function is continuous
          if (degree >= 1) {
            // test first derivative (central difference quotient)
            errorTest((basis.eval(l, i, x + dx) - basis.eval(l, i, x - dx)) /
                          (2.0 * dx),
                      basisEvalDx(basis, l, i, x), tol1);
          }

          if (degree >= 2) {
            // test second derivative (central difference quotient)
            errorTest((basisEvalDx(basis, l, i, x + dx) -
                       basisEvalDx(basis, l, i, x - dx)) /
                          (2.0 * dx),
                      basisEvalDxDx(basis, l, i, x), tol2);
          }
        }
      }

      BOOST_CHECK_LE(discontinuities, maxDiscontinuitiesCount);
    }
  }
}

void boundTest(SBasis& basis, level_t l, index_t i, double lowerBound,
               double upperBound) {
  for (size_t j = 0; j <= 100; j++) {
    const double x = static_cast<double>(j) / 100.0;
    const double fx = basis.eval(l, i, x);
    BOOST_CHECK_GE(fx, lowerBound);
    BOOST_CHECK_LE(fx, upperBound);
  }
}

void bsplinePropertiesTest(SBasis& basis, level_t start_level = 1,
                           bool modified = false) {
  // Test basic B-spline properties (mixed monotonicity, bounds) for
  // level >= start_level.
  const double tol = 1e-10;

  for (level_t l = 0; l < 6; l++) {
    const index_t hInv = static_cast<index_t>(1) << l;

    for (index_t i = 1; i < hInv; i += 2) {
      // test bounds
      const double upperBound =
          (((!modified) || ((i > 1) && (i < hInv - 1))) ? 1.0 : 2.02);
      boundTest(basis, l, i, -tol, upperBound);

      // rising at the beginning
      bool falling = false;
      double fx = std::numeric_limits<double>::quiet_NaN();

      for (size_t j = 0; j < 100; j++) {
        const double x = static_cast<double>(j) / 100.0;
        const double fxNew = basis.eval(l, i, x);

        if (!std::isnan(fx)) {
          if (falling) {
            // hope we're still falling
            BOOST_CHECK_LE(fxNew, fx + tol);
          } else if (fxNew < fx - tol) {
            // we're now falling (and weren't until now)
            falling = true;
          }
        }

        fx = fxNew;
      }
    }
  }
}

void fundamentalSplineTest(SBasis& basis, bool modified = false,
                           bool nak = false) {
  const level_t startLevel = 1;
  const double tol = 1e-10;

  for (level_t l = startLevel; l < 6; l++) {
    const index_t hInv = static_cast<index_t>(1) << l;

    for (index_t i = 1; i < hInv; i += 2) {
      // test bounds
      double lowerBound;
      double upperBound;

      if (modified) {
        lowerBound = -0.3;
        upperBound = (((i > 1) && (i < hInv - 1)) ? 1.0 : 2.3);
      } else if (nak) {
        lowerBound = -0.5;
        upperBound = 1.4;
      } else {
        lowerBound = -0.3;
        upperBound = 1.0;
      }

      boundTest(basis, l, i, lowerBound - tol, upperBound + tol);

      for (index_t i2 = 0; i2 <= hInv; i2++) {
        // test Lagrange property
        if ((!modified) || ((i > 1) && (i < hInv - 1)) ||
            ((i2 > 0) && (i2 < hInv))) {
          const double x = static_cast<double>(i2) / static_cast<double>(hInv);
          const double fx = basis.eval(l, i, x);

          BOOST_CHECK_SMALL(fx - ((i == i2) ? 1.0 : 0.0), 1e-10);
        }

        // test sign
        if (i2 < hInv) {
          const double sign =
              (((i2 - i) % 2 == 0) ? 1.0 : -1.0) * ((i2 < i) ? -1.0 : 1.0);

          for (size_t j = 0; j < 100; j++) {
            const double x =
                (static_cast<double>(i2) + static_cast<double>(j) / 100.0) /
                static_cast<double>(hInv);
            const double fx = basis.eval(l, i, x);

            if (sign == 1.0) {
              BOOST_CHECK_GE(fx, -1e-10);
            } else {
              BOOST_CHECK_LE(fx, 1e-10);
            }
          }
        }
      }
    }
  }
}

void weaklyFundamentalSplineTest(SBasis& basis, bool modified = false) {
  const level_t startLevel = 1;

  for (level_t l = startLevel; l < 6; l++) {
    const index_t hInv = static_cast<index_t>(1) << l;

    for (index_t i = 1; i < hInv; i += 2) {
      for (index_t i2 = 0; i2 <= hInv; i2 += 2) {
        // test Lagrange property
        if ((!modified) || ((i > 1) && (i < hInv - 1)) ||
            ((i2 > 0) && (i2 < hInv))) {
          const double x = static_cast<double>(i2) / static_cast<double>(hInv);
          const double fx = basis.eval(l, i, x);

          BOOST_CHECK_SMALL(fx - ((i == i2) ? 1.0 : 0.0), 1e-10);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE(TestAlgorithms)

BOOST_AUTO_TEST_CASE(TestLinearBasis) {
  sgpp::base::SLinearBase basis;
  linearUniformUnmodifiedTest(basis);
  derivativesTest(basis, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearBoundaryBasis) {
  sgpp::base::SLinearBoundaryBase basis;
  linearLevelZeroTest(basis);
  linearUniformUnmodifiedTest(basis);
  derivativesTest(basis, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearClenshawCurtisBasis) {
  sgpp::base::SLinearClenshawCurtisBoundaryBase basis;
  linearLevelZeroTest(basis);
  linearClenshawCurtisBoundaryTest(basis);
  derivativesTest(basis, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearModifiedBasis) {
  sgpp::base::SLinearModifiedBase basis;
  linearModifiedTest(basis);
  derivativesTest(basis, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearStretchedBasis) {
  sgpp::base::SLinearStretchedBase basis;

  const std::vector<double> p = {1.312631659, 0.96716821};
  const std::vector<double> pos0 = {0.5, 0.5};
  const std::vector<double> pos1 = {7.0, 7.0};
  const std::vector<double> testValues = {0.125020255230769, 0.071872032307692};

  stretchedBasisTest(basis, p, pos0, pos1, testValues);
}

BOOST_AUTO_TEST_CASE(TestBsplineBasis) {
  // Test B-spline Noboundary basis.
  sgpp::base::SBsplineBase basis(1);
  linearUniformUnmodifiedTest(basis);

  for (size_t p = 2; p <= 11; p++) {
    sgpp::base::SBsplineBase basis(p);
    bsplinePropertiesTest(basis);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineBoundaryBasis) {
  // Test B-spline Boundary basis.
  sgpp::base::SBsplineBoundaryBase basis(1);
  linearLevelZeroTest(basis);
  linearUniformUnmodifiedTest(basis);

  for (size_t p = 2; p <= 11; p++) {
    sgpp::base::SBsplineBoundaryBase basis(p);
    bsplinePropertiesTest(basis, 0);
    derivativesTest(basis, basis.getDegree() - 1, 0);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineClenshawCurtisBoundaryBasis) {
  // Test B-spline ClenshawCurtis basis.
  sgpp::base::SBsplineClenshawCurtisBase basis(1);
  linearLevelZeroTest(basis);
  linearClenshawCurtisBoundaryTest(basis);

  for (size_t p = 2; p <= 11; p++) {
    sgpp::base::SBsplineClenshawCurtisBase basis(p);
    bsplinePropertiesTest(basis, 0);
    derivativesTest(basis, basis.getDegree() - 1, 0);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineModifiedBasis) {
  // Test modified B-spline basis.
  sgpp::base::SBsplineModifiedBase basis(1);
  linearModifiedTest(basis);

  for (size_t p = 2; p <= 11; p++) {
    sgpp::base::SBsplineModifiedBase basis(p);
    bsplinePropertiesTest(basis, 1, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineClenshawCurtisModifiedBasis) {
  // Test modified B-spline ClenshawCurtis basis.
  for (size_t p = 2; p <= 11; p++) {
    sgpp::base::SBsplineModifiedClenshawCurtisBase basis(p);
    bsplinePropertiesTest(basis, 1, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestFundamentalNakSplineBasis) {
  // Test fundamental not-a-knot spline basis.
  const size_t pMax = 5;

  for (size_t p = 1; p <= pMax; p++) {
    sgpp::base::SFundamentalNakSplineBase basis(p);
    fundamentalSplineTest(basis, false, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestFundamentalSplineBasis) {
  // Test fundamental spline basis.
  const size_t pMax = 11;

  for (size_t p = 1; p <= pMax; p++) {
    sgpp::base::SFundamentalSplineBase basis(p);
    fundamentalSplineTest(basis);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestFundamentalSplineModifiedBasis) {
  // Test modified fundamental spline basis.

  // choose lower p_max, because 2nd derivative tests fail easily at high
  // degrees (finite differences not accurate enoguh)
  const size_t pMax = 7;

  for (size_t p = 1; p <= pMax; p++) {
    sgpp::base::SFundamentalSplineModifiedBase basis(p);
    fundamentalSplineTest(basis, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestWeaklyFundamentalNakSplineBasis) {
  // Test weakly fundamental not-a-knot spline basis.
  for (size_t p = 1; p <= 7; p++) {
    sgpp::base::SWeaklyFundamentalNakSplineBase basis(p);
    weaklyFundamentalSplineTest(basis);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestWeaklyFundamentalNakSplineModifiedBasis) {
  // Test modified weakly fundamental not-a-knot spline basis.
  for (size_t p = 1; p <= 7; p++) {
    sgpp::base::SWeaklyFundamentalNakSplineModifiedBase basis(p);
    weaklyFundamentalSplineTest(basis, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestWeaklyFundamentalSplineBasis) {
  // Test weakly fundamental spline basis.
  for (size_t p = 1; p <= 7; p++) {
    sgpp::base::SWeaklyFundamentalSplineBase basis(p);
    weaklyFundamentalSplineTest(basis);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestNakBsplineBasis) {
  // Test not-a-knot B-spline basis.
  for (size_t p = 1; p <= 7; p++) {
    sgpp::base::SNakBsplineBase basis(p);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestNakBsplineModifiedBasis) {
  // Test modified weakly fundamental not-a-knot spline basis.
  for (size_t p = 1; p <= 5; p++) {
    sgpp::base::SNakBsplineModifiedBase basis(p);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestWaveletBasis) {
  // Test Wavelet Noboundary basis.
  sgpp::base::SWaveletBase basis;
  derivativesTest(basis, 2, 1, 2);
}

BOOST_AUTO_TEST_CASE(TestWaveletBoundaryBasis) {
  // Test Wavelet Boundary basis.
  sgpp::base::SWaveletBoundaryBase basis;
  derivativesTest(basis, 2, 1, 2);
}

BOOST_AUTO_TEST_CASE(TestWaveletModifiedBasis) {
  // Test modified Wavelet basis.
  sgpp::base::SWaveletModifiedBase basis;
  derivativesTest(basis, 2, 1, 2);
}

BOOST_AUTO_TEST_CASE(TestGetAffectedBasisFunctions) {
  GridPoint i(1);
  GridStorage s(1);
  sgpp::base::SLinearBase b;

  i.set(0, 1, 1);
  s.insert(i);

  sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBase> ga(s);
  std::vector<std::pair<size_t, double>> x;
  DataVector y(1, 0.25);

  ga(b, y, x);

  BOOST_CHECK_EQUAL(x[0].first, 0U);
  BOOST_CHECK_EQUAL(x[0].second, 0.5);
}

BOOST_AUTO_TEST_CASE(TestGetAffectedBasisFunctionsBoundary) {
  GridPoint i(1);
  GridStorage s(1);
  sgpp::base::SLinearBoundaryBase b;

  i.set(0, 0, 0);
  s.insert(i);
  i.set(0, 0, 1);
  s.insert(i);
  i.set(0, 1, 1);
  s.insert(i);

  sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBoundaryBase> ga(s);
  std::vector<std::pair<size_t, double>> x;
  DataVector y(1, 0.5);

  ga(b, y, x);

  BOOST_CHECK_EQUAL(x[0].first, 0U);
  BOOST_CHECK_EQUAL(x[0].second, 0.5);

  BOOST_CHECK_EQUAL(x[1].first, 1U);
  BOOST_CHECK_EQUAL(x[1].second, 0.5);

  BOOST_CHECK_EQUAL(x[2].first, 2U);
  BOOST_CHECK_EQUAL(x[2].second, 1.0);
}

BOOST_AUTO_TEST_CASE(TestGetAffectedBasisFunctionsStretched) {
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;

  BoundingBox1D dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;
  Stretching stretch({dimBound}, {str1d});

  GridPoint i(1);
  GridStorage s(1);
  s.setStretching(stretch);

  sgpp::base::SLinearStretchedBoundaryBase b;

  i.set(0, 0, 0);
  s.insert(i);
  i.set(0, 0, 1);
  s.insert(i);
  i.set(0, 1, 1);
  s.insert(i);

  sgpp::base::GetAffectedBasisFunctions<
      sgpp::base::SLinearStretchedBoundaryBase> ga(s);
  std::vector<std::pair<size_t, double>> x;
  DataVector y(1, 0.25);

  ga(b, y, x);

  BOOST_CHECK_EQUAL(x[0].first, 0U);
  BOOST_CHECK_CLOSE(x[0].second, 1.0384615384615385, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()

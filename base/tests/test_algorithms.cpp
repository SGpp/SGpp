#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include "BasisEval.hpp"

using namespace SGPP;
using namespace SGPP::base;

const bool use_double_precision =
#if USE_DOUBLE_PRECISION
  true;
#else
  false;
#endif /* USE_DOUBLE_PRECISION */

void basisTest(SBasis& basis,
               const std::vector<level_t>& levels,
               const std::vector<index_t>& indices,
               const std::vector<SGPP::float_t>& points,
               const std::vector<SGPP::float_t>& testValues) {
  const size_t n = indices.size();
  BOOST_CHECK_EQUAL(levels.size(), n);
  BOOST_CHECK_EQUAL(points.size(), n);
  BOOST_CHECK_EQUAL(testValues.size(), n);

  for (size_t i = 0; i < n; i++) {
    SGPP::float_t val = basis.eval(levels[i], indices[i], points[i]);
    BOOST_CHECK_CLOSE(val, testValues[i], 1e-10);
  }
}

void stretchedBasisTest(SLinearStretchedBase& basis,
                        const std::vector<SGPP::float_t>& p,
                        const std::vector<SGPP::float_t>& pos0,
                        const std::vector<SGPP::float_t>& pos1,
                        const std::vector<SGPP::float_t>& testValues) {
  const size_t n = p.size();
  BOOST_CHECK_EQUAL(pos0.size(), n);
  BOOST_CHECK_EQUAL(pos1.size(), n);
  BOOST_CHECK_EQUAL(testValues.size(), n);

  for (size_t i = 0; i < n; i++) {
    SGPP::float_t val = basis.eval(p[i], pos0[i], pos1[i]);
    BOOST_CHECK_CLOSE(val, testValues[i], 1e-10);
  }
}

void linearUniformUnmodifiedTest(SBasis& basis) {
  const std::vector<level_t> levels = {1, 1, 1, 1, 1, 2, 2, 3, 3, 3};
  const std::vector<index_t> indices = {1, 1, 1, 1, 1, 1, 3, 1, 1, 1};

  const std::vector<SGPP::float_t> points = {0.5, 0.75, 0.875, 0.0, 1.0, 0.75,
                                             0.75, 0.0, 0.125, 0.25
                                            };

  const std::vector<SGPP::float_t> testValues = {1.0, 0.5, 0.25, 0.0, 0.0, 0.0,
                                                 1.0, 0.0, 1.0, 0.0
                                                };

  basisTest(basis, levels, indices, points, testValues);
}

void linearModifiedTest(SBasis& basis) {
  const std::vector<level_t> levels = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
                                       3, 3, 3
                                      };
  const std::vector<index_t> indices = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3,
                                        1, 1, 1
                                       };

  const std::vector<SGPP::float_t> points = {0.5, 0.75, 0.875, 0.0, 1.0, 0.0,
                                             0.125, 0.25, 0.375, 0.75, 0.75, 1.0, 0.0, 0.125, 0.25
                                            };

  const std::vector<SGPP::float_t> testValues = {1.0, 1.0, 1.0, 1.0, 1.0, 2.0,
                                                 1.5, 1.0, 0.5, 0.0, 1.0, 2.0, 2.0, 1.0, 0.0
                                                };

  basisTest(basis, levels, indices, points, testValues);
}

void linearLevelZeroTest(SBasis& basis) {
  // Test boundary linear basis functions (level 0).
  for (index_t i = 0; i <= 1; i++) {
    for (size_t j = 0; j <= 4; j++) {
      const SGPP::float_t x = static_cast<SGPP::float_t>(j) / 4.0;
      BOOST_CHECK_EQUAL(basis.eval(0, i, x), i * x + (1 - i) * (1.0 - x));
    }
  }
}

SGPP::float_t ccKnot(level_t l, index_t i) {
  // Return Clenshaw-Curtis knot with given level and index.
  return 0.5 * (std::cos(M_PI * (1.0 - static_cast<SGPP::float_t>(i) /
                                 std::pow(2.0, static_cast<SGPP::float_t>(l)))) + 1.0);
}

void linearClenshawCurtisTest(SBasis& basis) {
  const std::vector<level_t> levels = {1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3};
  const std::vector<level_t> indices = {1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 1, 1};

  const std::vector<SGPP::float_t> points = {
    0.5, 0.75, 0.875, 0.0, 1.0, 0.75,
    0.75, ccKnot(2, 3), 0.0, 0.125, ccKnot(3, 1), 0.25
  };

  const std::vector<double> testValuesDouble = {
      1.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.25 / (double(ccKnot(2, 3)) - 0.5), 1.0, 0.0,
      1.0 - (0.125 - double(ccKnot(3, 1))) / (double(ccKnot(3, 2) - ccKnot(3, 1))), 1.0, 0.0
    };

  const std::vector<SGPP::float_t> testValues( testValuesDouble.begin(),
					       testValuesDouble.end() );

  basisTest(basis, levels, indices, points, testValues);
}

void errorTest(SGPP::float_t x, SGPP::float_t y, SGPP::float_t tol) {
  if (std::abs(x) >= 10.0) {
    BOOST_CHECK_SMALL((x - y) / x, tol);
  } else {
    BOOST_CHECK_SMALL(x - y, SGPP::float_t(10.0) * tol);
  }
}

void derivativesTest(SBasis& basis, size_t degree = 2, level_t startLevel = 1,
                     size_t maxDiscontinuitiesCount = 0) {
  // Test derivatives (up to order deg, max. 2) of basis functions for
  // level >= start_level. Allow for max_discontinuities_count
  // discontinuities (e.g. for Wavelets which are cut off).

  if (!use_double_precision) {
    // skip test when using single precision, because then
    // the derivatives are not exact enough
    return;
  }

  const SGPP::float_t dx = 1e-8;
  const SGPP::float_t tol1 = 1e-3;
  const SGPP::float_t tol2 = 1e-2;
  const SGPP::float_t discontinuityTol = 1e5;

  // levels
  for (level_t l = 0; l < 6; l++) {
    // indices
    for (index_t i = 1; i < (static_cast<index_t>(1) << l); i += 2) {
      if (degree >= 1) {
        // test first derivative at boundary (central difference quotient)
        errorTest((basis.eval(l, i, 2.0 * dx) - basis.eval(l, i, 0.0)) /
                  (2.0 * dx), basisEvalDx(basis, l, i, dx), tol1);
        errorTest((basis.eval(l, i, 1.0) - basis.eval(l, i, 1.0 - 2.0 * dx)) /
                  (2.0 * dx), basisEvalDx(basis, l, i, 1.0 - dx), tol1);
      }

      if (degree >= 2) {
        // test second derivative at boundary (central difference quotient)
        errorTest((basisEvalDx(basis, l, i, 2.0 * dx) -
                   basisEvalDx(basis, l, i, 0.0)) /
                  (2.0 * dx), basisEvalDxDx(basis, l, i, dx), tol2);
        errorTest((basisEvalDx(basis, l, i, 1.0) -
                   basisEvalDx(basis, l, i, 1.0 - 2.0 * dx)) /
                  (2.0 * dx), basisEvalDxDx(basis, l, i, 1.0 - dx), tol2);
      }

      size_t discontinuities = 0;

      for (size_t j = 1; j < 100; j++) {
        const SGPP::float_t x = static_cast<SGPP::float_t>(j) / 100.0;

        if (std::abs(basis.eval(l, i, x + dx) - basis.eval(l, i, x - dx)) >
            discontinuityTol * dx) {
          // discontinuity found
          discontinuities++;
        } else {
          // test derivatives only if the function is continuous
          if (degree >= 1) {
            // test first derivative (central difference quotient)
            errorTest((basis.eval(l, i, x + dx) - basis.eval(l, i, x - dx)) /
                      (2.0 * dx), basisEvalDx(basis, l, i, x), tol1);
          }

          if (degree >= 2) {
            // test second derivative (central difference quotient)
            errorTest((basisEvalDx(basis, l, i, x + dx) -
                       basisEvalDx(basis, l, i, x - dx)) /
                      (2.0 * dx), basisEvalDxDx(basis, l, i, x), tol2);
          }
        }
      }

      BOOST_CHECK_LE(discontinuities, maxDiscontinuitiesCount);
    }
  }
}

void boundTest(SBasis& basis, level_t l, index_t i,
               SGPP::float_t lowerBound, SGPP::float_t upperBound) {
  for (size_t j = 0; j <= 100; j++) {
    const SGPP::float_t x = static_cast<SGPP::float_t>(j) / 100.0;
    const SGPP::float_t fx = basis.eval(l, i, x);
    BOOST_CHECK_GE(fx, lowerBound);
    BOOST_CHECK_LE(fx, upperBound);
  }
}

void bsplinePropertiesTest(SBasis& basis, level_t start_level = 1,
                           bool modified = false) {
  // Test basic B-spline properties (mixed monotonicity, bounds) for
  // level >= start_level.
  const SGPP::float_t tol =
#if USE_DOUBLE_PRECISION
    0.0;
#else
    1e-4;
#endif

  for (level_t l = 0; l < 6; l++) {
    const index_t hInv = static_cast<index_t>(1) << l;

    for (index_t i = 1; i < hInv; i += 2) {
      // test bounds
      const SGPP::float_t upperBound =
        (((!modified) || ((i > 1) && (i < hInv - 1))) ? 1.0 : 2.02);
      boundTest(basis, l, i, 0.0, upperBound);

      // rising at the beginning
      bool falling = false;
      SGPP::float_t fx = NAN;

      for (size_t j = 0; j < 100; j++) {
        const SGPP::float_t x = static_cast<SGPP::float_t>(j) / 100.0;
        const SGPP::float_t fxNew = basis.eval(l, i, x);

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

void fundamentalSplineTest(SBasis& basis, bool modified = false) {
  const level_t startLevel = 1;

  for (level_t l = startLevel; l < 6; l++) {
    const index_t hInv = static_cast<index_t>(1) << l;

    for (index_t i = 1; i < hInv; i += 2) {
      // test bounds
      const SGPP::float_t upperBound =
        (((!modified) || ((i > 1) && (i < hInv - 1))) ? 1.0 : 2.3);
      boundTest(basis, l, i, -0.3, upperBound + 1e-10);

      for (index_t i2 = 0; i2 <= hInv; i2++) {
        // test Lagrange property
        if ((!modified) || ((i > 1) && (i < hInv - 1)) ||
            ((i2 > 0) && (i2 < hInv))) {
          const SGPP::float_t x = static_cast<SGPP::float_t>(i2) /
                                  static_cast<SGPP::float_t>(hInv);
          const SGPP::float_t fx = basis.eval(l, i, x);
#if USE_DOUBLE_PRECISION
          BOOST_CHECK_SMALL(fx - ((i == i2) ? 1.0 : 0.0), 1e-10);
#else
          BOOST_CHECK_SMALL(fx - ((i == i2) ? 1.0 : 0.0), 1e-2);
#endif
        }

        // test sign
        if (i2 < hInv) {
          const SGPP::float_t sign = (((i2 - i) % 2 == 0) ? 1.0 : -1.0) *
                                     ((i2 < i) ? -1.0 : 1.0);

          for (size_t j = 0; j < 100; j++) {
            const SGPP::float_t x = (static_cast<SGPP::float_t>(i2) +
                                     static_cast<SGPP::float_t>(j) / 100.0) /
                                    static_cast<SGPP::float_t>(hInv);
            const SGPP::float_t fx = basis.eval(l, i, x);

            if (sign == 1.0) {
#if USE_DOUBLE_PRECISION
              BOOST_CHECK_GE(fx, -1e-10);
#else
              BOOST_CHECK_GE(fx, -1e-10);
#endif
            } else {
#if USE_DOUBLE_PRECISION
              BOOST_CHECK_LE(fx, 1e-10);
#else
              BOOST_CHECK_LE(fx, 1e-10);
#endif
            }
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestLinearBasis) {
  SLinearBase basis;
  linearUniformUnmodifiedTest(basis);
  derivativesTest(basis, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearBoundaryBasis) {
  SLinearBoundaryBase basis;
  linearLevelZeroTest(basis);
  linearUniformUnmodifiedTest(basis);
  derivativesTest(basis, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearClenshawCurtisBasis) {
  SLinearClenshawCurtisBase basis;
  linearLevelZeroTest(basis);
  linearClenshawCurtisTest(basis);
  derivativesTest(basis, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearModifiedBasis) {
  SLinearModifiedBase basis;
  linearModifiedTest(basis);
  derivativesTest(basis, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestLinearStretchedBasis) {
  SLinearStretchedBase basis;

  const std::vector<SGPP::float_t> p = {1.312631659, 0.96716821};
  const std::vector<SGPP::float_t> pos0 = {0.5, 0.5};
  const std::vector<SGPP::float_t> pos1 = {7.0, 7.0};
  const std::vector<SGPP::float_t> testValues = {0.125020255230769,
                                                 0.071872032307692
                                                };

  stretchedBasisTest(basis, p, pos0, pos1, testValues);
}

BOOST_AUTO_TEST_CASE(TestBsplineBasis) {
  // Test B-spline Noboundary basis.
  SBsplineBase basis(1);
  linearUniformUnmodifiedTest(basis);

  const size_t pMax = (use_double_precision ? 11 : 6);

  for (size_t p = 2; p <= pMax; p++) {
    SBsplineBase basis(p);
    bsplinePropertiesTest(basis);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineBoundaryBasis) {
  // Test B-spline Boundary basis.
  SBsplineBoundaryBase basis(1);
  linearLevelZeroTest(basis);
  linearUniformUnmodifiedTest(basis);

  const size_t pMax = (use_double_precision ? 11 : 6);

  for (size_t p = 2; p <= pMax; p++) {
    SBsplineBoundaryBase basis(p);
    bsplinePropertiesTest(basis, 0);
    derivativesTest(basis, basis.getDegree() - 1, 0);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineClenshawCurtisBasis) {
  // Test B-spline ClenshawCurtis basis.
  SBsplineClenshawCurtisBase basis(1);
  linearLevelZeroTest(basis);
  linearClenshawCurtisTest(basis);

  const size_t pMax = (use_double_precision ? 11 : 6);

  for (size_t p = 2; p <= pMax; p++) {
    SBsplineClenshawCurtisBase basis(p);
    bsplinePropertiesTest(basis, 0);
    derivativesTest(basis, basis.getDegree() - 1, 0);
  }
}

BOOST_AUTO_TEST_CASE(TestBsplineModifiedBasis) {
  // Test modified B-spline basis.
  SBsplineModifiedBase basis(1);
  linearModifiedTest(basis);

  const size_t pMax = (use_double_precision ? 11 : 6);

  for (size_t p = 2; p <= pMax; p++) {
    SBsplineModifiedBase basis(p);
    bsplinePropertiesTest(basis, 1, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestFundamentalSplineBasis) {
  // Test fundamental spline basis.
  const size_t pMax = 11;

  for (size_t p = 1; p <= pMax; p++) {
    SFundamentalSplineBase basis(p);
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
    SFundamentalSplineModifiedBase basis(p);
    fundamentalSplineTest(basis, true);
    derivativesTest(basis, basis.getDegree() - 1);
  }
}

BOOST_AUTO_TEST_CASE(TestWaveletBasis) {
  // Test Wavelet Noboundary basis.
  SWaveletBase basis;
  derivativesTest(basis, 2, 1, 2);
}

BOOST_AUTO_TEST_CASE(TestWaveletBoundaryBasis) {
  // Test Wavelet Boundary basis.
  SWaveletBoundaryBase basis;
  derivativesTest(basis, 2, 1, 2);
}

BOOST_AUTO_TEST_CASE(TestWaveletModifiedBasis) {
  // Test modified Wavelet basis.
  SWaveletModifiedBase basis;
  derivativesTest(basis, 2, 1, 2);
}

BOOST_AUTO_TEST_CASE(TestGetAffectedBasisFunctions) {
  GridIndex i(1);
  GridStorage s(1);
  SLinearBase b;

  i.set(0, 1, 1);
  s.insert(i);

  GetAffectedBasisFunctions<SLinearBase> ga(&s);
  std::vector<std::pair<size_t, SGPP::float_t>> x;
  DataVector y(1, 0.25);

  ga(b, y, x);

  BOOST_CHECK_EQUAL(x[0].first, 0);
  BOOST_CHECK_EQUAL(x[0].second, 0.5);
}

BOOST_AUTO_TEST_CASE(TestGetAffectedBasisFunctionsBoundary) {
  GridIndex i(1);
  GridStorage s(1);
  SLinearBoundaryBase b;

  i.set(0, 0, 0);
  s.insert(i);
  i.set(0, 0, 1);
  s.insert(i);
  i.set(0, 1, 1);
  s.insert(i);

  GetAffectedBasisFunctions<SLinearBoundaryBase> ga(&s);
  std::vector<std::pair<size_t, SGPP::float_t>> x;
  DataVector y(1, 0.5);

  ga(b, y, x);

  BOOST_CHECK_EQUAL(x[0].first, 0);
  BOOST_CHECK_EQUAL(x[0].second, 0.5);

  BOOST_CHECK_EQUAL(x[1].first, 1);
  BOOST_CHECK_EQUAL(x[1].second, 0.5);

  BOOST_CHECK_EQUAL(x[2].first, 2);
  BOOST_CHECK_EQUAL(x[2].second, 1.0);
}

BOOST_AUTO_TEST_CASE(TestGetAffectedBasisFunctionsStretched) {
  Stretching1D str1d;
  str1d.type = "log";
  str1d.x_0 = 1;
  str1d.xsi = 10;

  DimensionBoundary dimBound;
  dimBound.leftBoundary = 0.5;
  dimBound.rightBoundary = 7;
  Stretching stretch(1, &dimBound, &str1d);

  GridIndex i(1);
  GridStorage s(1);
  s.setStretching(stretch);

  SLinearStretchedBoundaryBase b;

  i.set(0, 0, 0);
  s.insert(i);
  i.set(0, 0, 1);
  s.insert(i);
  i.set(0, 1, 1);
  s.insert(i);

  GetAffectedBasisFunctions<SLinearStretchedBoundaryBase> ga(&s);
  std::vector<std::pair<size_t, SGPP::float_t>> x;
  DataVector y(1, 0.25);

  ga(b, y, x);

  BOOST_CHECK_EQUAL(x[0].first, 0);
  BOOST_CHECK_CLOSE(x[0].second, 1.0384615384615385, 1e-5);
}

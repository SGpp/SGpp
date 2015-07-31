#define BOOST_TEST_DYN_LINK
#include <vector>
#include <boost/test/unit_test.hpp>

#include "sgpp/base/operation/hash/common/basis/LinearBasis.hpp"
#include "sgpp/base/operation/hash/common/basis/Basis.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/grid/GridStorage.hpp"

using namespace SGPP::base;

struct F {
  F() {};
  ~F() {};

  void baseTest( SBasis& basis,
                 std::vector< level_t >& levels,
                 std::vector< index_t >& indices,
                 std::vector< SGPP::float_t >& points,
                 std::vector< SGPP::float_t >& testvals ) {
    BOOST_CHECK( levels.size() == indices.size() );
    BOOST_CHECK( points.size() == indices.size() );
    BOOST_CHECK( testvals.size() == indices.size() );

    for ( size_t i = 0; i < levels.size(); ++i ) {
      SGPP::float_t val = basis.eval( levels[i], indices[i], points[i] );

      // here it is possible to set different tolerances for float and double
#if USE_DOUBLE_PRECISION
      BOOST_CHECK_CLOSE( val, testvals[i], 0.0 );
#else
      BOOST_CHECK_CLOSE( val, testvals[i], 0.0 );
#endif
    }
  }

  void derivativesTest( SBasis& b,
                        unsigned int deg = 2,
                        level_t start_level = 1,
                        unsigned int max_discontinuities_count = 0 ) {
    /*

    # levels
    for l in range(start_level, 6):
      # indices
      for i in range(1, 2**l, 2):
        f = lambda y: b.eval(l, i, y)

        if deg >= 1:
          # test first derivative at boundary (central difference quotient)
          df = lambda y: b.evalDx(l, i, y)
          self.errorTest((f(2.0*dx) - f(0.0)) / (2.0*dx), df(dx), tol1)
          self.errorTest((f(1.0) - f(1.0-2.0*dx)) / (2.0*dx), df(1.0-dx), tol1)
          #self.assertAlmostEqual((f(2.0*dx) - f(0.0)) / (2.0*dx), df(dx), delta=tol1)
          #self.assertAlmostEqual((f(1.0) - f(1.0-2.0*dx)) / (2.0*dx), df(1.0-dx), delta=tol1)
        if deg >= 2:
          # test second derivative at boundary (central difference quotient)
          ddf = lambda y: b.evalDxDx(l, i, y)
          self.errorTest((df(2.0*dx) - df(0.0)) / (2.0*dx), ddf(dx), tol2)
          self.errorTest((df(1.0) - df(1.0-2.0*dx)) / (2.0*dx), ddf(1.0-dx), tol2)
          #self.assertAlmostEqual((df(2.0*dx) - df(0.0)) / (2.0*dx), ddf(dx), delta=tol2)
          #self.assertAlmostEqual((df(1.0) - df(1.0-2.0*dx)) / (2.0*dx), ddf(1.0-dx), delta=tol2)

        # count discontinuities
        discontinuities = 0
        for x in [j/100.0 for j in range(1, 100)]:
          if abs(f(x+dx) - f(x-dx)) > discontinuityTol * dx:
            # discontinuity found
            discontinuities += 1
          else:
            # test derivatives only if the function is continuous
            if deg >= 1:
              # test first derivative (central difference quotient)
              self.errorTest((f(x+dx) - f(x-dx)) / (2.0*dx), df(x), tol1)
              #self.assertAlmostEqual((f(x+dx) - f(x-dx)) / (2.0*dx), df(x), delta=tol1)
            if deg >= 2:
              # test second derivative (central difference quotient)
              self.errorTest((df(x+dx) - df(x-dx)) / (2.0*dx), ddf(x), tol2)
              #self.assertAlmostEqual((df(x+dx) - df(x-dx)) / (2.0*dx), ddf(x), delta=tol2)
        self.assertLessEqual(discontinuities, max_discontinuities_count)
     */

    /* Test derivatives (up to order deg, max. 2) of basis functions for
     * level >= start_level. Allow for max_discontinuities_count discontinuities
     * (e.g. for Wavelets which are cut off).
     */

    // skip test when using single precision, because then
    // the derivatives are not exact enough
#ifndef USING_DOUBLE_PRECISION
    return;
#endif

    SGPP::float_t dx = 1e-8;
    SGPP::float_t tol1 = 1e-3;
    SGPP::float_t tol2 = 1e-2;
    SGPP::float_t discontinuityTol = 1e5;

    // f = lambda y: b.eval(l, i, y)
    // df = lambda y: b.evalDx(l, i, y)

    // levels
    for ( level_t l = start_level; l < 6; ++l ) {
      // indices
      for ( level_t i = 1; i < static_cast<level_t>( std::pow(2, l) ); i += 2 ) {
        if ( deg >= 1 ) {
          //test first derivative at boundary (central difference quotient)
          /*
          errorTest(
            ( b.eval( l, i, 2.0*dx ) - b.eval( l, i, 0.0) ) / (2.0*dx),
            b.evalDx( l, i, dx ), tol1 );*/
          //errorTest( (f(1.0) - f(1.0-2.0*dx)) / (2.0*dx), df(1.0-dx), tol1 );
        }



      }

    }
  }

  void errorTest( SGPP::float_t x, SGPP::float_t y, SGPP::float_t tol ) {
    if ( std::abs(x) >= 10.0 ) {
      BOOST_CHECK( std::abs( x - y) / std::abs(x) < tol );
    } else {
      BOOST_CHECK_SMALL( std::abs(x - y), 10.0 * tol );
    }
  }

  void linearUniformUnmodifiedTest( SBasis& basis ) {
    std::vector< level_t > levels = { 1, 1, 1, 1, 1, 2, 2, 3, 3, 3 };
    std::vector< index_t > indices = { 1, 1, 1, 1, 1, 1, 3, 1, 1, 1 };

    std::vector< SGPP::float_t > points = { 0.5, 0.75, 0.875, 0.0, 1.0,
                                            0.75, 0.75, 0.0, 0.125, 0.25
                                          };

    std::vector< SGPP::float_t > testvals { 1.0, 0.5, 0.25, 0.0, 0.0,
                                            0.0, 1.0, 0.0, 1.0, 0.0 };

    baseTest( basis, levels, indices, points, testvals );
  }
};

BOOST_AUTO_TEST_SUITE(Algorithms)

BOOST_FIXTURE_TEST_SUITE( TestBase, F )

BOOST_AUTO_TEST_CASE(testLinear) {
  /*
   *   def testLinear(self):
    from pysgpp import SLinearBase
    b = SLinearBase()
    self.linearUniformUnmodifiedTest(b)
    self.derivativesTest(b, deg=0)
   */
  SLinearBase b;

  linearUniformUnmodifiedTest( b );
}
// end test suite testbase
BOOST_AUTO_TEST_SUITE_END()

// end test suite algorithms
BOOST_AUTO_TEST_SUITE_END()



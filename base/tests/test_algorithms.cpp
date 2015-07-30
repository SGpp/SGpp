#define BOOST_TEST_DYN_LINK
#include <vector>
#include <boost/test/unit_test.hpp>

#include "sgpp/base/operation/hash/common/basis/LinearBasis.hpp"
#include "sgpp/base/operation/hash/common/basis/Basis.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"

using namespace SGPP::base;

struct F {
  F(){};
  ~F(){};

  /*
   *   def baseTest(self, b, points):
  from pysgpp import cvar
  for t in points:
    val = b.eval(t[0], t[1], t[2])
    errorMessage = "%f != %f => (%d, %d) @ %f"%(val, t[3], t[0], t[1], t[2])
    if cvar.USING_DOUBLE_PRECISION:
      self.assertAlmostEqual(val, t[3], msg=errorMessage)
    else:
      self.assertAlmostEqual(val, t[3], msg=errorMessage, places=6)
   */
  void baseTest( SBasis& basis,
		 std::vector< unsigned int >& levels,
		 std::vector< unsigned int >& indices,
		 std::vector< SGPP::float_t >& points,
		 std::vector< SGPP::float_t >& testvals ){
    BOOST_CHECK( levels.size() == indices.size() );
    BOOST_CHECK( points.size() == indices.size() );
    BOOST_CHECK( testvals.size() == indices.size() );

    for( size_t i=0; i<levels.size(); ++i ){
      SGPP::float_t val = basis.eval( levels[i], indices[i], points[i] );

      // here it is possible to set different tolerances for float and double
#if USE_DOUBLE_PRECISION
      BOOST_CHECK_CLOSE( val, testvals[i], 0.0 );
#else
      BOOST_CHECK_CLOSE( val, testvals[i], 0.0 );
#endif
    }
  }

  void linearUniformUnmodifiedTest( SBasis& basis ){
    std::vector< unsigned int > levels = { 1, 1, 1, 1, 1, 2, 2, 3, 3, 3 };
    std::vector< unsigned int > indices = { 1, 1, 1, 1, 1, 1, 3, 1, 1, 1 };

    std::vector< SGPP::float_t > points = { 0.5, 0.75, 0.875, 0.0, 1.0,
					    0.75, 0.75, 0.0, 0.125, 0.25 };

    std::vector< SGPP::float_t > testvals { 1.0, 0.5, 0.25, 0.0, 0.0,
					    0.0, 1.0, 0.0, 1.0, 0.0 };

    baseTest( basis, levels, indices, points, testvals );
  }
};

BOOST_AUTO_TEST_SUITE(Algorithms)

BOOST_FIXTURE_TEST_SUITE( TestBase, F )

BOOST_AUTO_TEST_CASE(testLinear){
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



#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


//declare testcase with name "Simple"
BOOST_AUTO_TEST_CASE(Simple) {

  double var1 = 7.0;
  double var2 = 8.0;

  BOOST_CHECK(var1 == var2);

}

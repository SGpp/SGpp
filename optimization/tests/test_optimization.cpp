#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BoostTestOptimization
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#include <omp.h>

struct GlobalFixture {
  const int MAX_NUM_THREADS = 2;

  GlobalFixture() {
    // limit number of OpenMP threads, because if chosen too high,
    // then testing will be very slow on machines with a high CPU load
    // (e.g., bad for Jenkins)
    if (omp_get_max_threads() > MAX_NUM_THREADS) {
      omp_set_num_threads(MAX_NUM_THREADS);
    }
  }
};

#if BOOST_VERSION >= 105900
BOOST_GLOBAL_FIXTURE(GlobalFixture);
#else
BOOST_GLOBAL_FIXTURE(GlobalFixture)
#endif /* BOOST_VERSION >= 105900 */

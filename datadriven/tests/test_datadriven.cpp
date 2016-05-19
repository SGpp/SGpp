// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DatadrivenTestModule
#include <boost/test/unit_test.hpp>
#include <omp.h>

// struct GlobalFixture {
//  const int MAX_NUM_THREADS = 2;
//
//  GlobalFixture() {
//    // limit number of OpenMP threads, because if chosen too high,
//    // then testing will be very slow on machines with a high CPU load
//    // (e.g., bad for Jenkins)
//    if (omp_get_max_threads() > MAX_NUM_THREADS) {
//      omp_set_num_threads(MAX_NUM_THREADS);
//    }
//  }
//};
//
// BOOST_GLOBAL_FIXTURE(GlobalFixture)

// fix for clang (from https://stackoverflow.com/a/33755176)
#ifdef __clang__
#include <string>

namespace boost {
namespace unit_test {
namespace ut_detail {

std::string normalize_test_case_name(const_string name) {
    return ((name[0] == '&') ? std::string(name.begin() + 1, name.size() - 1) :
                               std::string(name.begin(), name.size()));
}

}  // namespace ut_detail
}  // namespace unit_test
}  // namespace boost
#endif

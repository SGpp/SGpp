// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppCombigridModule

#include <boost/test/unit_test.hpp>

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

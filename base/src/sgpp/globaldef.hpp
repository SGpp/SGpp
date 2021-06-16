// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_GLOBALDEF_HPP_
#define SGPP_GLOBALDEF_HPP_

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

/* Define math constants.
 * These seem not to be part of the C++ standard,
 * and they're missing, e.g., on MinGW-w64.
 */
#ifdef _MSC_VER
#undef _MATH_DEFINES_DEFINED
#define _MATH_DEFINES_DEFINED
#endif
// e
#undef M_E
#define M_E 2.7182818284590452354
// log_10(2)
#undef M_LOG2E
#define M_LOG2E 1.4426950408889634074
// log_10(e)
#undef M_LOG10E
#define M_LOG10E 0.43429448190325182765
// ln(2)
#undef M_LN2
#define M_LN2 0.69314718055994530942
// ln(10)
#undef M_LN10
#define M_LN10 2.30258509299404568402
// pi
#undef M_PI
#define M_PI 3.14159265358979323846
// pi/2
#undef M_PI_2
#define M_PI_2 1.57079632679489661923
// pi/4
#undef M_PI_4
#define M_PI_4 0.78539816339744830962
// 1/pi
#undef M_1_PI
#define M_1_PI 0.31830988618379067154
// 2/pi
#undef M_2_PI
#define M_2_PI 0.63661977236758134308
// 2/sqrt(pi)
#undef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
// sqrt(2)
#undef M_SQRT2
#define M_SQRT2 1.41421356237309504880
// 1/sqrt(2)
#undef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440

// 1/sqrt(2 pi)
#undef M_1_SQRT2PI
#define M_1_SQRT2PI 0.398942280401432702863218082712

#if __cplusplus == 201103L
#include <memory>  // NOLINT(build/include)

namespace std {
// Implementation for "make_unique" in c++11 as
// it doesn't contain this function.  (see std::make_shared)
// This function is part of the C++14 (and newer) standard.
template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {  // NOLINT(build/c++11)
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
}
#endif

#if __cplusplus == 201402L
// Hotfix for c++14, which seems to have stricter requirements for the availability of the size_t
// type.
#include <cstddef>
#include <memory>  // NOLINT(build/include)
#endif

#endif /* SGPP_GLOBALDEF_HPP_ */

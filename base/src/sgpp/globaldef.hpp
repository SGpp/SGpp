// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_GLOBALDEF_HPP_
#define SGPP_GLOBALDEF_HPP_

////////////////////////////////////////////////////////////////////////////////
///////////         GLOBAL DEFINES OF THE SG++ PROJECT             /////////////
////////////////////////////////////////////////////////////////////////////////
/// The library can be built to support double or single precision floating  ///
/// point numbers, but only one at a time.                                   ///
/// In order to use both in the same implementation, two separate libraries  ///
/// need to be built and linked.                                             ///
/// The separation of the symbols is achieved through the use of different   ///
/// top level namespaces: sg for the double precision library, sgsp for the  ///
/// single precision library.                                                ///
/// This is implemented via a preprocessor switch that changes the top level ///
/// namespace and defines "USE_DOUBLE_PRECISION" to be used for distinction  ///
/// between single and double precision implementations of vectorized        ///
/// kernels.                                                                 ///
////////////////////////////////////////////////////////////////////////////////
///////////         PLEASE READ AGAIN AND UNDERSTAND!              /////////////
////////////////////////////////////////////////////////////////////////////////

#ifndef USE_DOUBLE_PRECISION
#define USE_DOUBLE_PRECISION 1
#endif

#if USE_DOUBLE_PRECISION==1
#define SGPP sg
#else
#define SGPP sgsp
#endif

/* Define math constants.
 * These seem not to be part of the C++ standard,
 * and they're missing, e.g., on MinGW-w64.
 */
#ifdef _MSC_VER
#undef  _MATH_DEFINES_DEFINED
#define _MATH_DEFINES_DEFINED
#endif
// e
#undef  M_E
#define M_E         2.7182818284590452354
// log_10(2)
#undef  M_LOG2E
#define M_LOG2E     1.4426950408889634074
// log_10(e)
#undef  M_LOG10E
#define M_LOG10E    0.43429448190325182765
// ln(2)
#undef  M_LN2
#define M_LN2       0.69314718055994530942
// ln(10)
#undef  M_LN10
#define M_LN10      2.30258509299404568402
// pi
#undef  M_PI
#define M_PI        3.14159265358979323846
// pi/2
#undef  M_PI_2
#define M_PI_2      1.57079632679489661923
// pi/4
#undef  M_PI_4
#define M_PI_4      0.78539816339744830962
// 1/pi
#undef  M_1_PI
#define M_1_PI      0.31830988618379067154
// 2/pi
#undef  M_2_PI
#define M_2_PI      0.63661977236758134308
// 2/sqrt(pi)
#undef  M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257390
// sqrt(2)
#undef  M_SQRT2
#define M_SQRT2     1.41421356237309504880
// 1/sqrt(2)
#undef  M_SQRT1_2
#define M_SQRT1_2   0.70710678118654752440

namespace SGPP {

#if USE_DOUBLE_PRECISION
/// the library uses 64 bit floating point numbers
typedef double float_t;
#else
/// the library uses 32 bit floating point numbers
typedef float float_t;
#endif

}
/* namespace SGPP */

#if __cplusplus == 201103L
#include <memory>

namespace std {
//Implementation for "make_unique" in c++11 as it doesn't contain this function.  (see std::make_shared)
//This function is part of the C++14 (and newer) standard.
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique(Args&& ...args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}
#endif

#endif /* SGPP_GLOBALDEF_HPP_ */

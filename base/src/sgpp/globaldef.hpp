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

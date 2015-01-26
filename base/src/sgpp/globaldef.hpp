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

#ifndef SGPP

#define USE_DOUBLE_PRECISION 1
#define SGPP sg

#else

#define USE_DOUBLE_PRECISION 0
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

#endif /* SGPP_GLOBALDEF_HPP_ */
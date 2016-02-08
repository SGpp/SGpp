/* ****************************************************************************
 * Copyright (C) 2010 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Janos Benk
#ifndef COMBIGRID_ULTILS_HPP_
#define COMBIGRID_ULTILS_HPP_

/** In this hpp file we define the utilities for the combination technique */
#include <math.h>
#include <iostream>
#include <typeinfo>
#include <stack>
#include <functional>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <complex>
#include <assert.h>
#include <float.h>


#define ABS(x) x < 0 ? -x : x
#define UNDEFINED_OPERATION_EXCEPTION -1
#define UNSUPPORTED_OPERATION_EXCEPTION -2
#define VIOLATED_PRECEDENCE_RELATION_EXCEPTION -3

/**
 * define macros for positive and negative infinities
 */
#define p_INF DBL_MAX
#define n_INF -(DBL_MAX)

// include open mp for shared memory programming
#include <omp.h>

// useful for checking the fullgrid vs coefficient template type precedence.. .

/**
 * Utility function that checks if the precedence of template type _FGTp is lower than that of
 * template type _CFTp. the function execetutes successfully if
 *
 * _CFTp >= _FGTp (i.e. _CFTp has higher precedence than _FGTp)
 *
 *
 * Use Case 1: check_type_precedence<float,double>(); --> should pass without exception (as _FGTp <= _CFTp)
 * Use Case 2: check_type_precedence<float,float>(); --> should pass without exception (as _FGTp = _CFTp)
 * Use Case 3: check_type_precedence<double,float>(); --> should throw an exception (as _FGTp > _CFTp )
 *
 */
using namespace std;

template<typename _FGTp, typename _CFTp>
void check_type_precedence() {

  if (typeid(_FGTp) == typeid(_CFTp)) {
    return;
  }

  if (typeid(_FGTp).before(typeid(_CFTp))) {
    throw VIOLATED_PRECEDENCE_RELATION_EXCEPTION;
  }

  /// which is not supposed to happen
}






//vector declarations only for internal use , to make the interface more simpler
typedef std::vector<double> DVector;

typedef std::vector<int> IVector;

#define COMBIGRID_OUT_WRN(str,FILE,LINE){\
  std::cout<<"\nWARNING: (in FILE:" << FILE << "; LINE:" << LINE << ") :\n" << str << "\n";}

#define COMBIGRID_OUT_ERR(str,FILE,LINE){\
  std::cout<<"FILE:" << FILE << " LINE:" << LINE << ": \n ERROR:" << str << "\n";}

#define COMBIGRID_OUT(str) { \
      std::cout << str << std::endl; };

#define COMBIGRID_ERROR_TEST(test,str) { \
      if (!(test)) std::cout << std::endl << "ERROR: " << str << std::endl; \
      assert(test); } \

#define COMBIGRID_ERROR_TEST_EQUAL(val1,val2,tolerance,str) { \
      if (fabs(val1-val2) > tolerance ) {std::cout << std::endl << "ERROR: " << str << \
                   " , v1:" << val1 << " , v2:" << val2 << " , tol:" << tolerance <<std::endl; \
      assert(false); }} \

#define COMBIGRID_ERROR_EXIT(str) { \
      std::cout << std::endl << "ERROR: " << str << std::endl; \
      assert(false); } \

/** write error message out*/
#define COMBIGRID_ERROR_MSG(str) { \
      std::cout << std::endl << "ERROR: " << str << std::endl; } \

/** Conditionally write out */
#define COMBIGRID_OUT_LEVEL(limit,level,str) \
     if (limit < level ) COMBIGRID_OUT(str);

/** Conditionally write out */
#define COMBIGRID_OUT_LEVEL1(level,str) \
     if (level > 1) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL2(level,str) \
     if (level > 2) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL3(level,str) \
     if (level > 3) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL4(level,str) \
     if (level > 4) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL5(level,str) \
     if (level > 5) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL6(level,str) \
     if (level > 6) COMBIGRID_OUT(str);

#define COMBIGRID_OUT_LEVEL7(level,str) \
     if (level > 7) COMBIGRID_OUT(str);

namespace combigrid {

/** vector with power two */
const int powerOfTwo[30] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
                             4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576,
                             2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
                             268435456
                           };

/** vector with one over power two */
const double oneOverPowOfTwo[30] = { 1.0 / 1.0, 1.0 / 2.0, 1.0 / 4.0, 1.0 / 8.0,
                                     1.0 / 16.0, 1.0 / 32.0, 1.0 / 64.0, 1.0 / 128.0, 1.0 / 256.0, 1.0
                                     / 512.0, 1.0 / 1024.0, 1.0 / 2048.0, 1.0 / 4096.0, 1.0 / 8192.0,
                                     1.0 / 16384.0, 1.0 / 32768.0, 1.0 / 65536.0, 1.0 / 131072.0, 1.0
                                     / 262144.0, 1.0 / 524288.0, 1.0 / 1048576.0, 1.0 / 2097152.0,
                                     1.0 / 4194304.0, 1.0 / 8388608.0, 1.0 / 16777216.0, 1.0 / 33554432.0,
                                     1.0 / 67108864.0, 1.0 / 134217728.0, 1.0 / 268435456.0
                                   };

/** the maximum double value */
inline double COMBIGRID_DMAX(double v1, double v2) {
  return (v1 < v2) ? v2 : v1;
}

/** the maximum int value */
inline int COMBIGRID_IMAX(int v1, int v2) {
  return (v1 < v2) ? v2 : v1;
}

/** the minimum double value */
inline double COMBIGRID_DMIN(double v1, double v2) {
  return (v1 > v2) ? v2 : v1;
}

/** the minimum int value */
inline int COMBIGRID_IMIN(int v1, int v2) {
  return (v1 > v2) ? v2 : v1;
}

/** the function C_{N}^K , combination of N,K*/
static int combination(int n, int k) {
  if ((k == 0) || (n == k))
    return 1;
  else if ((k == 1) || (n == k + 1))
    return n;
  else if ((k == 2) || (n == k + 2))
    return n * (n - 1) / 2;
  else
    return combination(n - 1, k) + combination(n - 1, k - 1);
}

/** calculates the L2 norm of one vector*/
inline double l2_norm(std::vector<double>* v1) {
  double diff = 0.0;

  for (unsigned int i = 0; i < v1->size(); i++) {
    diff = diff + v1->at(i) * v1->at(i);
  }

  return sqrt(diff / double(v1->size()));
}

/** calculates the Inf norm of one vector*/
inline double inf_norm(std::vector<double>* v1) {
  double diff = 0.0, tmp;

  for (unsigned int i = 0; i < v1->size(); i++) {
    tmp = fabs(v1->at(i));
    diff = (tmp > diff) ? tmp : diff;
  }

  return diff;
}

/** multiply two vectors, result will be in the first vector */
inline void vect_mul(std::vector<double>* v1, std::vector<double>* v2) {
  COMBIGRID_ERROR_TEST( v1->size() == v2->size(),
                        " vect_mul , size do not match v1->size():" << v1->size() << " , v2->size():" <<
                        v2->size());

  for (unsigned int i = 0; i < v1->size(); i++) {
    v1->at(i) = v1->at(i) * v2->at(i);
  }
}

/** difference of two vectors, result will be in the first vector */
inline void vect_diff(std::vector<double>* v1, std::vector<double>* v2) {
  COMBIGRID_ERROR_TEST( v1->size() == v2->size(),
                        " vect_diff , size do not match v1->size():" << v1->size() << " , v2->size():"
                        << v2->size());

  for (unsigned int i = 0; i < v1->size(); i++) {
    v1->at(i) = v1->at(i) - v2->at(i);
  }
}

/** v1 = coefv1*v1 + coefv2 * v2 */
inline void vect_add_mul(double coefv1, std::vector<double>* v1, double coefv2,
                         std::vector<double>* v2) {
  COMBIGRID_ERROR_TEST( v1->size() == v2->size(),
                        " vect_add_mul , size do not match v1->size():" << v1->size() <<
                        " , v2->size():" << v2->size());

  for (unsigned int i = 0; i < v1->size(); i++) {
    v1->at(i) = coefv1 * v1->at(i) + coefv2 * v2->at(i);
  }
}

/** computes the scalar product of two vectors */
inline void scalar_product(std::vector<double>* v1, std::vector<double>* v2,
                           double& result) {
  COMBIGRID_ERROR_TEST( v1->size() == v2->size(),
                        " vect_add_mul , size do not match v1->size():" << v1->size() <<
                        " , v2->size():" << v2->size());
  result = 0.0;

  for (unsigned int i = 0; i < v1->size(); i++) {
    result = result + v1->at(i) * v2->at(i);
  }
}

/** sets the values of the vector to a given value */
inline void vect_setvalue(std::vector<double>* v1, double newValue) {
  for (unsigned int i = 0; i < v1->size(); i++) {
    v1->at(i) = newValue;
  }
}

/** plot one vector */
inline void plot_vect(int level, int verb, std::vector<double>* v1,
                      const char* stri) {
  if (verb > level) {
    std::cout << stri << "=[" << v1->at(0);

    for (unsigned int i = 1; i < v1->size(); i++) {
      std::cout << "," << v1->at(i);
    }

    std::cout << "];" << std::endl;
  }
}
}

#endif /* COMBIGRID_ULTILS_HPP_ */

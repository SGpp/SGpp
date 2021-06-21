// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BASIS_HPP
#define BASIS_HPP

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Basis class for basis functions.
 */
template <class LT, class IT>
class Basis {
 public:
  /**
   * Evaluate the basis function with given level and index.
   *
   * @param level   level of the basis function
   * @param index   index of the basis function
   * @param x       evaluation point
   * @result        value of the basis function.
   */
  virtual double eval(LT level, IT index, double x) = 0;

  /**
   * Evaluate the basis functions derivative with given level and index.
   *
   * @param level   level of the basis function
   * @param index   index of the basis function
   * @param x       evaluation point
   * @result        value of the basis functions derivative.
   */
  virtual double evalDx(LT level, IT index, double x) = 0;

  /**
   * Returns the polynomial degree of the basis
   *
   * @return polynomial degree of the basis
   */
  virtual size_t getDegree() const = 0;

  /**
   * returns the integal of the current basis function
   *
   * @param level   level of the basis function
   * @param index   index of the basis function
   * @return
   */
  virtual double getIntegral(LT level, IT index) = 0;

  /**
   * Destructor.
   */
  virtual ~Basis() {}

  // Basis();
  /*private:
  Basis(Basis const&);
  Basis& operator=(Basis const&);
  */
};

typedef Basis<unsigned int, unsigned int> SBasis;

}  // namespace base
}  // namespace sgpp

#endif  // BASIS_HPP

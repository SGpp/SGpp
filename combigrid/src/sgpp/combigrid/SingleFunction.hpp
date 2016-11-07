// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_

#include <sgpp/globaldef.hpp>

#include <functional>

namespace sgpp {
namespace combigrid {

/**
 * Wrapper class for std::function<double(double)>. This is necessary because SWIG can't
 * handle std::function objects properly.
 */
class SingleFunction {
 public:
  typedef std::function<double(double)> function_type;

 private:
  function_type func;

 public:
  /**
   * for function pointers
   */
  explicit SingleFunction(double (*ptr)(double));

  /**
   * for lambdas or function objects
   */
  template <typename T>
  explicit SingleFunction(T f) : func(f) {}

  /**
   * Default constructor, creating a constant zero function.
   */
  SingleFunction() : func([](double x) { return 0.0; }) {}

  /**
   * Evaluates the function.
   */
  double operator()(double param);

  /**
   * Evaluates the function (for python use etc., does the same as operator()).
   */
  double call(double param);

  /**
   * Returns the corresponding std::function<double(double)> object.
   */
  function_type getLambdaExpression();
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_ */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_MULTIFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_MULTIFUNCTION_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <functional>

namespace sgpp {
namespace combigrid {

/**
 * Wrapper for std::function<double(base::DataVector const &)>. This is necessary because SWIG can't
 * handle std::function objects properly.
 */
class MultiFunction {
 public:
  typedef std::function<double(base::DataVector const &)> function_type;

 private:
  function_type func;

 public:
  /**
   * for function pointers
   */
  explicit MultiFunction(double (*ptr)(base::DataVector const &));

  /**
   * for lambdas or function objects
   */
  template <typename T>
  explicit MultiFunction(T const &f) : func(f) {}

  /**
   * Evaluates the function.
   */
  double operator()(base::DataVector const &vec) const;

  /**
   * Evaluates the function (for python use etc., does the same as operator()).
   */
  double call(base::DataVector const &vec) const;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_MULTIFUNCTION_HPP_ */

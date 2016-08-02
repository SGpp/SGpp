/*
 * SingleFunction.hpp
 *
 *  Created on: 27.02.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_

#include <sgpp/globaldef.hpp>
#include <functional>

namespace sgpp {
namespace combigrid {

class SingleFunction {
 public:
  typedef std::function<double(double)> function_type;

 private:
  function_type func;

 public:
  /**
   * for function pointers
   */
  SingleFunction(double (*ptr)(double));

  /**
   * for lambdas or function objects
   */
  template <typename T>
  explicit SingleFunction(T f)
      : func(f) {}

  double operator()(double param);
  double call(double param);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_ */

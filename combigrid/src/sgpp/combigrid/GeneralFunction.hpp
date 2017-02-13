// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GENERALFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GENERALFUNCTION_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <functional>

namespace sgpp {
namespace combigrid {

/**
 * Wrapper for std::function<Out(In)>. This is necessary because SWIG can't
 * handle std::function objects properly.
 */
template <typename Out, typename In>
class GeneralFunction {
 public:
  typedef std::function<Out(In)> function_type;
  // typedef In input_type;
  typedef Out output_type;

 private:
  function_type func;

 public:
  /**
   * for function pointers
   */
  explicit GeneralFunction(Out (*ptr)(In)) : func(ptr) {}

  /**
   * for lambdas or function objects
   */
  template <typename T>
  explicit GeneralFunction(T const &f) : func(f) {}

  /**
   * Default constructor, creating a constant zero function.
   */
  GeneralFunction() : func([](In) { return Out(); }) {}

  /**
   * Evaluates the function.
   */
  Out operator()(In in) const { return func(in); }

  /**
   * Evaluates the function (for python use etc., does the same as operator()).
   */
  Out call(In in) const { return func(in); }

  function_type getStdFunction() const { return func; }
};

template <typename Out>
class GeneralFunction1 {
 public:
  typedef std::function<Out()> function_type;
  // typedef In input_type;
  typedef Out output_type;

 private:
  function_type func;

 public:
  /**
   * for function pointers
   */
  explicit GeneralFunction1(Out (*ptr)()) : func(ptr) {}

  /**
   * for lambdas or function objects
   */
  template <typename T>
  explicit GeneralFunction1(T const &f) : func(f) {}

  /**
   * Default constructor, creating a constant zero function.
   */
  GeneralFunction1() : func([]() { return Out(); }) {}

  /**
   * Evaluates the function.
   */
  Out operator()() const { return func(); }

  /**
   * Evaluates the function (for python use etc., does the same as operator()).
   */
  Out call() const { return func(); }

  function_type getStdFunction() const { return func; }
};

typedef GeneralFunction<double, base::DataVector const &> MultiFunction;
typedef GeneralFunction<double, double> SingleFunction;

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GENERALFUNCTION_HPP_ */

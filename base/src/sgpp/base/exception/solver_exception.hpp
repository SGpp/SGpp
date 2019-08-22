// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SOLVER_EXCEPTION_HPP
#define SOLVER_EXCEPTION_HPP

#include <sgpp/globaldef.hpp>

#include <exception>
#include <cstddef>


namespace sgpp {
namespace base {

/**
 * Exception that is thrown in case of a solver operation failure
 *
 */
class solver_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit solver_exception(const char* msg) noexcept : msg(msg) {
  }

  /**
   * Standard Constructor
   */
  solver_exception() noexcept : msg(nullptr) { }

  /**
   * Destructor
   */
  ~solver_exception() noexcept override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  const char* what() const noexcept override {
    if (msg) {
      return msg;
    } else {
      return "solver_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;
};

}  // namespace base
}  // namespace sgpp

#endif /* SOLVER_EXCEPTION_HPP */

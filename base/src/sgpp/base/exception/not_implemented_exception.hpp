// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <exception>

namespace sgpp {
namespace base {

/**
 * Exception that is thrown in case of an application failure
 *
 */
class not_implemented_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit not_implemented_exception(const char* msg) noexcept : msg(msg) {}

  /**
   * Standard Constructor
   */
  not_implemented_exception() noexcept : msg(NULL) {}

  /**
   * Destructor
   */
  ~not_implemented_exception() noexcept override {}

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  const char* what() const noexcept override {
    if (msg) {
      return msg;
    } else {
      return "not_implemented_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;
};

}  // namespace base
}  // namespace sgpp

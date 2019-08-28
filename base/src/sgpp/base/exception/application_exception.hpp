// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef APPLICATION_EXCEPTION_HPP
#define APPLICATION_EXCEPTION_HPP

#include <sgpp/globaldef.hpp>

#include <exception>
#include <cstddef>


namespace sgpp {
namespace base {

/**
 * Exception that is thrown in case of an application failure
 *
 */
class application_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit application_exception(const char* msg) noexcept : msg(msg) {
  }

  /**
   * Standard Constructor
   */
  application_exception() noexcept : msg(nullptr) { }

  /**
   * Destructor
   */
  ~application_exception() noexcept override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  const char* what() const noexcept override {
    if (msg) {
      return msg;
    } else {
      return "application_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;
};

}  // namespace base
}  // namespace sgpp

#endif /* APPLICATION_EXCEPTION_HPP */

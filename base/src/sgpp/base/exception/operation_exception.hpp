// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATION_EXCEPTION_HPP
#define OPERATION_EXCEPTION_HPP

#include <sgpp/globaldef.hpp>

#include <exception>
#include <cstddef>
#include <string>


namespace sgpp {
namespace base {

/**
 * Exception that is thrown in case of a grid operation failure
 *
 */
class operation_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit operation_exception(const char* msg) noexcept : msg(msg) {
  }

  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit operation_exception(std::string msg) noexcept : stringMsg(msg),
    msg(stringMsg.c_str()) {
  }

  /**
   * Standard Constructor
   */
  operation_exception() noexcept : msg(nullptr) { }

  /**
   * Destructor
   */
  ~operation_exception() noexcept override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  const char* what() const noexcept override {
    if (msg) {
      return msg;
    } else {
      return "operation_exception: general failure";
    }
  }

 protected:
  /// the exception message
  std::string stringMsg;
  /// the exception message as C string
  const char* msg;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATION_EXCEPTION_HPP */

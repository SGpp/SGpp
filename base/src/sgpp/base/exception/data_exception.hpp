// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATA_EXCEPTION_HPP
#define DATA_EXCEPTION_HPP

#include <sgpp/globaldef.hpp>

#include <exception>
#include <cstddef>


namespace sgpp {
namespace base {

/**
 * Exception that is thrown in case of a data failure (conversion, creation, ...).
 *
 */
class data_exception : public std::exception {
 public:
  /**
   * Create a new exception (constructor) with some message.
   *
   * @param msg The exception message
   */
  explicit data_exception(const char* msg) noexcept : msg(msg) {
  }

  /**
   * Create default exception (constructor).
   */
  data_exception() noexcept : msg(NULL) { }

  /**
   * Destructor
   */
  ~data_exception() noexcept override { }

  /**
   * Return message of exception object.
   *
   * @return Returns the message specified in the constructor, otherwise a general text
   */
  const char* what() const noexcept override {
    if (msg) {
      return msg;
    } else {
      return "data_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;
};

}  // namespace base
}  // namespace sgpp

#endif /* DATA_EXCEPTION_HPP */

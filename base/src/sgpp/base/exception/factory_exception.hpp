// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FACTORY_EXCEPTION_HPP
#define FACTORY_EXCEPTION_HPP

#include <sgpp/globaldef.hpp>

#include <exception>
#include <cstddef>


namespace sgpp {
namespace base {

/**
 * Exception that is thrown in case of a grid failure
 *
 */
class factory_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit factory_exception(const char* msg) noexcept : msg(msg) {
  }

  /**
   * Standard Constructor
   */
  factory_exception() noexcept : msg(NULL) { }

  /**
   * Destructor
   */
  ~factory_exception() noexcept override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  const char* what() const noexcept override {
    if (msg) {
      return msg;
    } else {
      return "factory_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;
};

}  // namespace base
}  // namespace sgpp

#endif /* FACTORY_EXCEPTION_HPP */

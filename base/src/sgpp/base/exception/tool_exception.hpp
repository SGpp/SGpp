// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TOOL_EXCEPTION_HPP
#define TOOL_EXCEPTION_HPP

#include <sgpp/globaldef.hpp>

#include <exception>
#include <cstddef>


namespace SGPP {
namespace base {

/**
 * Exception that is thrown in case of a tool operation failure
 *
 */
class tool_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  explicit tool_exception(const char* msg) throw() : msg(msg) {
  }

  /**
   * Standard Constructor
   */
  tool_exception() throw() : msg(NULL) { }

  /**
   * Destructor
   */
  ~tool_exception() throw() override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  const char* what() const throw() override {
    if (msg) {
      return msg;
    } else {
      return "tool_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;
};

}  // namespace base
}  // namespace SGPP

#endif /* TOOL_EXCEPTION_HPP */

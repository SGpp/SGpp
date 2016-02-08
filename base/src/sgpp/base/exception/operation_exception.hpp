// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATION_EXCEPTION_HPP
#define OPERATION_EXCEPTION_HPP

#include <exception>
#include <cstddef>
#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
  operation_exception(const char* msg) throw() : msg(msg) {
  }

  /**
   * Constructor
   *
   * @param msg the exception message
   */
  operation_exception(std::string msg) throw() : stringMsg(msg),
    msg(stringMsg.c_str()) {
  }

  /**
   * Standard Constructor
   */
  operation_exception() throw() : msg(NULL) { }

  /**
   * Destructor
   */
  virtual ~operation_exception() throw() override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  virtual const char* what() const throw() override {
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
}  // namespace SGPP

#endif /* OPERATION_EXCEPTION_HPP */

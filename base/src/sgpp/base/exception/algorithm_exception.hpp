// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORITHM_EXCEPTION_HPP
#define ALGORITHM_EXCEPTION_HPP

#include <exception>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Exception that is thrown in case of an application failure
 *
 */
class algorithm_exception : public std::exception {
 public:
  /**
   * Constructor
   *
   * @param msg the exception message
   */
  algorithm_exception(const char* msg) throw() : msg(msg) {
  }

  /**
   * Standard Constructor
   */
  algorithm_exception() throw() : msg(NULL) { }

  /**
   * Destructor
   */
  virtual ~algorithm_exception() throw() override { }

  /**
   * throw method that have to be implemented
   *
   * @return returns the message specified in the constructor otherwise a general text
   */
  virtual const char* what() const throw() override {
    if (msg) {
      return msg;
    } else {
      return "algorithm_exception: general failure";
    }
  }

 protected:
  /// the exception message
  const char* msg;

};

}  // namespace base
}  // namespace SGPP

#endif /* ALGORITHM_EXCEPTION_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FILE_EXCEPTION_HPP
#define FILE_EXCEPTION_HPP

#include <exception>
#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Exception that is thrown in case of a file failure
     *
     * @version $HEAD$
     */
    class file_exception : public std::exception {
      public:
        /**
         * Constructor
         *
         * @param msg the exception message
         */
        file_exception(const char* msg) throw() : msg(msg) {
        }

        /**
         * Standard Constructor
         */
        file_exception() throw() : msg(NULL) { }

        /**
         * Destructor
         */
        virtual ~file_exception() throw() { }

        /**
         * throw method that have to be implemented
         *
         * @return returns the message specified in the constructor otherwise a general text
         */
        virtual const char* what() const throw() {
          if (msg) {
            return msg;
          } else {
            return "file_exception: general failure";
          }
        }

      protected:
        /// the exception message
        const char* msg;
    };

  }
}

#endif /* FILE_EXCEPTION_HPP */
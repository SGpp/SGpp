// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef APPLICATION_EXCEPTION_HPP
#define APPLICATION_EXCEPTION_HPP

#include <exception>
#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Exception that is thrown in case of an application failure
     *
     * @version $HEAD$
     */
    class application_exception : public std::exception {
      public:
        /**
         * Constructor
         *
         * @param msg the exception message
         */
        application_exception(const char* msg) throw() : msg(msg) {
        }

        /**
         * Standard Constructor
         */
        application_exception() throw() : msg(NULL) { }

        /**
         * Destructor
         */
        virtual ~application_exception() throw() { }

        /**
         * throw method that have to be implemented
         *
         * @return returns the message specified in the constructor otherwise a general text
         */
        virtual const char* what() const throw() {
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

  }
}

#endif /* APPLICATION_EXCEPTION_HPP */
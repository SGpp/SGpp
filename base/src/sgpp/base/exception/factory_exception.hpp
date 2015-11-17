// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef FACTORY_EXCEPTION_HPP
#define FACTORY_EXCEPTION_HPP

#include <exception>
#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        factory_exception(const char* msg) throw() : msg(msg) {
        }

        /**
         * Standard Constructor
         */
        factory_exception() throw() : msg(NULL) { }

        /**
         * Destructor
         */
        virtual ~factory_exception() throw() override { }

        /**
         * throw method that have to be implemented
         *
         * @return returns the message specified in the constructor otherwise a general text
         */
        virtual const char* what() const throw() override {
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

  }
}

#endif /* FACTORY_EXCEPTION_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef GENERATION_EXCEPTION_HPP
#define GENERATION_EXCEPTION_HPP

#include <exception>
#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Exception that is thrown in case of a grid generation failure
     *
     * @version $HEAD$
     */
    class generation_exception : public std::exception {
      public:
        /**
         * Constructor
         *
         * @param msg the exception message
         */
        generation_exception(const char* msg) throw() : msg(msg) {
        }

        /**
         * Standared Constructor
         */
        generation_exception() throw() : msg(NULL) { }

        /**
         * Destructor
         */
        virtual ~generation_exception() throw() { }

        /**
         * throw method that have to be implemented
         *
         * @return returns the message specified in the constructor otherwise a general text
         */
        virtual const char* what() const throw() {
          if (msg) {
            return msg;
          } else {
            return "generation_exception: failure generating grid";
          }
        }

      protected:
        /// the exception message
        const char* msg;
    };

  }
}

#endif /* GENERATION_EXCEPTION_HPP */
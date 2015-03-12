// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALPARTIALDERIVATIVE_HPP
#define OPERATIONNAIVEEVALPARTIALDERIVATIVE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace base {

    /**
     * Abstract operation for evaluating a partial derivative of a linear combination of
     * basis functions.
     * The "naive" is indicating that classes implementing this operation should use a "naive"
     * approach, e.g. by evaluating all basis functions by brute force.
     */
    class OperationNaiveEvalPartialDerivative {
      public:
        /**
         * Constructor.
         */
        OperationNaiveEvalPartialDerivative() {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalPartialDerivative() {
        }

        /**
         * Pure virtual method for evaluating a partial derivative of a linear combination of
         * basis functions.
         *
         * @param alpha     coefficient vector
         * @param point     evaluation point
         * @param derivDim  dimension in which the partial derivative should be taken (0, ..., d-1)
         * @return          value of the partial derivative of the linear combination
         */
        float_t evalPartialDerivative(DataVector& alpha,
                                      std::vector<float_t>& point,
                                      size_t derivDim) {
          DataVector p(point);
          return evalPartialDerivative(alpha, p, derivDim);
        }

        /**
         * Convenience function for using DataVector as points.
         *
         * @param alpha     coefficient vector
         * @param point     evaluation point
         * @param derivDim  dimension in which the partial derivative should be taken (0, ..., d-1)
         * @return          value of the partial derivative of the linear combination
         */
        virtual float_t evalPartialDerivative(DataVector& alpha,
                                              DataVector& point,
                                              size_t derivDim) = 0;
    };

  }
}

#endif /* OPERATIONNAIVEEVALPARTIALDERIVATIVE_HPP */

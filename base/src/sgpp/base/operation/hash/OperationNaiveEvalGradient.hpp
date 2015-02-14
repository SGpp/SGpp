// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALGRADIENT_HPP
#define OPERATIONNAIVEEVALGRADIENT_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
  namespace base {

    /**
     * Abstract operation for evaluating a linear combination of basis functions and its gradient.
     * The "naive" is indicating that classes implementing this operation should use a "naive"
     * approach, e.g. by evaluating all basis functions by brute force.
     */
    class OperationNaiveEvalGradient {
      public:
        /**
         * Constructor.
         */
        OperationNaiveEvalGradient() {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalGradient() {
        }

        /**
         * Pure virtual method for evaluating a linear combination of basis functions and its gradient.
         *
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient vector of the linear combination
         * @return                  value of the linear combination
         */
        virtual float_t evalGradient(base::DataVector& alpha, const std::vector<float_t>& point,
                                    base::DataVector& gradient) = 0;

        /**
         * Convenience function for using base::DataVector as points.
         *
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient vector of the linear combination
         * @return                  value of the linear combination
         */
        virtual float_t evalGradient(base::DataVector& alpha, base::DataVector& point,
                                    base::DataVector& gradient) {
          const std::vector<float_t> p(point.getPointer(), point.getPointer() + point.getSize());
          return evalGradient(alpha, p, gradient);
        }
    };

  }
}

#endif /* OPERATIONNAIVEEVALGRADIENT_HPP */

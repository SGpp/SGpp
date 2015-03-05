// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALHESSIAN_HPP
#define OPERATIONNAIVEEVALHESSIAN_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace SGPP {
  namespace base {

    /**
     * Abstract operation for evaluating a linear combination of basis functions, its gradient
     * and its Hessian.
     * The "naive" is indicating that classes implementing this operation should use a "naive"
     * approach, e.g. by evaluating all basis functions by brute force.
     */
    class OperationNaiveEvalHessian {
      public:
        /**
         * Constructor.
         */
        OperationNaiveEvalHessian() {
        }

        /**
         * Virtual destructor.
         */
        virtual ~OperationNaiveEvalHessian() {
        }

        /**
         * Pure virtual method for evaluating a linear combination of basis functions, its gradient
         * and its Hessian.
         *
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient vector of the linear combination
         * @param[out]  hessian     Hessian matrix of the linear combination
         * @return                  value of the linear combination
         */
        virtual float_t evalHessian(base::DataVector& alpha, const std::vector<float_t>& point,
                                    base::DataVector& gradient, base::DataMatrix& hessian) = 0;

        /**
         * Convenience function for using base::DataVector as points.
         *
         * @param       alpha       coefficient vector
         * @param       point       evaluation point
         * @param[out]  gradient    gradient vector of the linear combination
         * @param[out]  hessian     Hessian matrix of the linear combination
         * @return                  value of the linear combination
         */
        virtual float_t evalHessian(base::DataVector& alpha, base::DataVector& point,
                                    base::DataVector& gradient, base::DataMatrix& hessian) {
          const std::vector<float_t> p(point.getPointer(), point.getPointer() + point.getSize());
          return evalHessian(alpha, p, gradient, hessian);
        }
    };

  }
}

#endif /* OPERATIONNAIVEEVALHESSIAN_HPP */

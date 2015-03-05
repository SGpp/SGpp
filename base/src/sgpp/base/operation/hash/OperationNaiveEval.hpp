// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVAL_HPP
#define OPERATIONNAIVEEVAL_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Operation that evaluates the current sparse grid function defined
     * by the coefficient vector @em alpha at a given point.
     * In contrast to OperationEval, implementations of OperationNaiveEval should
     * use a "naive" method for evaluating sparse grid functions, e.g. evaluate
     * all basis functions by brute force.
     *
     * @todo (pflueged) Use eval(DataVector& alpha, DataVector& point) as default
     */
    class OperationNaiveEval {
      public:
        /**
         * Default constructor.
         */
        OperationNaiveEval() {}

        /**
         * Destructor
         */
        virtual ~OperationNaiveEval() {}

        /**
         * Evaluates the sparse grid function at a given point.
         *
         * @param alpha The coefficients of the sparse grid's basis functions
         * @param point The coordinates of the evaluation point
         */
        virtual float_t eval(DataVector& alpha, std::vector<float_t>& point) = 0;

        /**
         * Evaluates the sparse grid function at a given point.
         *
         * @param alpha The coefficients of the sparse grid's basis functions
         * @param point The coordinates of the evaluation point
         */
        float_t eval(DataVector& alpha, DataVector& point) {
          std::vector<float_t> p;

          for (size_t i = 0; i < point.getSize(); i++) {
            p.push_back(point[i]);
          }

          return eval(alpha, p);
        }
    };

  }
}

#endif /* OPERATIONNAIVEEVAL_HPP */

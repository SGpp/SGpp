// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATURELINEAR_HPP
#define OPERATIONQUADRATURELINEAR_HPP

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Quadrature on sparse grid, linear grid without boundaries
     */
    class OperationQuadratureLinear : public OperationQuadrature {
      public:
        /**
         * Constructor of OperationQuadratureLinear
         *
         * @param storage Pointer to the grid's GridStorage object
         */
        OperationQuadratureLinear(GridStorage* storage) : storage(storage) {}

        virtual ~OperationQuadratureLinear() {}

        /**
         * Quadrature for piecewise linear hat basis functions. Computes
         * @f[ \sum_{\vec{l}} 2^{-|\vec{l}|}\alpha_{\vec{l}}. @f]
         *
         * @param alpha Coefficient vector for current grid
         */
        virtual double doQuadrature(DataVector& alpha);

      protected:
        // Pointer to the grid's GridStorage object
        GridStorage* storage;
    };

  }
}

#endif /* OPERATIONQUADRATURE_HPP */
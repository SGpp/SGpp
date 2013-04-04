/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONQUADRATURELINEAR_HPP
#define OPERATIONQUADRATURELINEAR_HPP

#include "base/operation/OperationQuadrature.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
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

/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)

#ifndef OPERATIONQUADRATURELINEARBOUND_HPP
#define OPERATIONQUADRATURELINEARBOUND_HPP

#include "base/operation/OperationQuadrature.hpp"
#include "base/grid/Grid.hpp"

namespace sg {
  namespace base {

    /**
     * Quadrature on sparse grid, linear grid without boundaries
     */
    class OperationQuadratureLinearBoundary : public OperationQuadrature {
      public:
        /**
         * Constructor of OperationQuadratureLinear
         *
         * @param storage Pointer to the grid's GridStorage object
         */
        OperationQuadratureLinearBoundary(GridStorage* storage) : storage(storage) {}

        virtual ~OperationQuadratureLinearBoundary() {}

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

#endif /* OPERATIONQUADRATURELINEARBOUND_HPP */

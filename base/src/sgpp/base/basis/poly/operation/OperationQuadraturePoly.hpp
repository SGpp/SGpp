// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREPOLY_HPP
#define OPERATIONQUADRATUREPOLY_HPP

#include <sgpp/base/operation/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/basis/poly/PolyBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Quadrature on sparse grid, polynomial grid without boundaries
     */
    class OperationQuadraturePoly : public OperationQuadrature {
      public:
        /**
         * Constructor of OperationQuadraturePoly
         *
         * @param storage Pointer to the grid's GridStorage object
         * @param degree the polynom's max. degree
         */
        OperationQuadraturePoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

        virtual ~OperationQuadraturePoly() {}

        /**
         * Quadrature for piecewise polynomial basis functions of max. degree 3
         *
         * @param alpha Coefficient vector for current grid
         */
        virtual double doQuadrature(DataVector& alpha);

      protected:
        // Pointer to the grid's GridStorage object
        GridStorage* storage;
        /// Poly Basis object
        SPolyBase base;

    };

  }
}

#endif /* OPERATIONQUADRATURE_HPP */
/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de)
// @author Benjamin

#ifndef OPERATIONFIRSTMOMENTLINEAR_HPP
#define OPERATIONFIRSTMOMENTLINEAR_HPP

#include <sgpp/base/operation/OperationFirstMoment.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * FirstMomemnt of sparse grid function, linear grid without boundaries
     */
    class OperationFirstMomentLinear : public OperationFirstMoment {
      public:
        /**
         * Constructor of OperationFirstMomentLinear
         *
         * @param storage Pointer to the grid's GridStorage object
         */
        OperationFirstMomentLinear(GridStorage* storage) : storage(storage) {}

        virtual ~OperationFirstMomentLinear() {}

        /**
         * Compute first moment of the function
         * @f[ \int_{\Omega} x\cdot f(x) dx. @f]
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

#endif /* OPERATIONFIRSTMOMENTLINEAR_HPP */

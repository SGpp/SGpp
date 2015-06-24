// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONLAPLACELINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implementation of Laplace for linear functions with boundaries
     *
     */
    class OperationLaplaceLinearStretchedBoundary: public UpDownOneOpDim {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's SGPP::base::GridStorage object
         */
        OperationLaplaceLinearStretchedBoundary(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceLinearStretchedBoundary();


      protected:
        virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACELINEARBOUNDARYSTRETCHED_HPP */
// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLAPLACELINEARSTRETCHED_HPP
#define OPERATIONLAPLACELINEARSTRETCHED_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    /**
     * Implementation for linear functions of Laplace Operation, linear grids without boundaries
     *
     * @version $HEAD$
     */
    class OperationLaplaceLinearStretched: public UpDownOneOpDim {
      public:
        /**
         * Constructor of OperationLaplaceLinearStretched
         *
         * @param storage Pointer to the grid's gridstorage obejct
         */
        OperationLaplaceLinearStretched(SGPP::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLaplaceLinearStretched();

      protected:
        virtual void specialOP(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim, size_t gradient_dim);

        virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

        virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACELINEARSTRETCHED_HPP */
// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSUBSPACES_HPP
#define DMSYSTEMMATRIXSUBSPACES_HPP

#include <string>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

//#include "AbstractOperationMultipleEval.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Class that implements the virtual class SGPP::base::OperationMatrix for the
     * application of classification for the Systemmatrix
     *
     * The Identity matrix is used as regularization operator.
     *
     * For the Operation B's mult and mutlTransposed functions
     * vectorized formulations are used.
     */
    class SystemMatrixLeastSquaresIdentity: public SGPP::datadriven::DMSystemMatrixBase {
      private:
        /// vectorization mode
        //ComputeKernelType kernelType;
        /// Number of orignal training instances
        size_t instances;
        /// Number of patched and used training instances
        size_t paddedInstances;
        /// OperationB for calculating the data matrix
        //AbstractOperationMultipleEval* B;
        SGPP::base::OperationMultipleEval* B;

        SGPP::base::Grid& grid;

        SGPP::datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to SGPP::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         */
        SystemMatrixLeastSquaresIdentity(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrix& trainData, float_t lambda);

        /**
         * Std-Destructor
         */
        virtual ~SystemMatrixLeastSquaresIdentity();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

        virtual void generateb(SGPP::base::DataVector& classes, SGPP::base::DataVector& b);

        virtual void prepareGrid();

        void setImplementation(SGPP::datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
          this->implementationConfiguration = operationConfiguration;
          this->B = SGPP::op_factory::createOperationMultipleEval(this->grid, *(this->dataset_), this->implementationConfiguration);
        }
    };

  }
}

#endif /* DMSYSTEMMATRIXSUBSPACES_HPP */

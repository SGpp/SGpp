// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITYMPI_HPP
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITYMPI_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/datadriven/algorithm/DMSystemMatrixBaseSP.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <string>

#include <sgpp/globaldef.hpp>

#if USE_DOUBLE_PRECISION==0


namespace SGPP {
  namespace parallel {

    /**
     * Class that implements the virtual class SGPP::base::OperationMatrix for the
     * application of classification for the Systemmatrix
     *
     * The Identity matrix is used as regularization operator.
     *
     * For the Operation B's mult and mutlTransposed functions
     * vectorized formulations are used.
     */
    class DMSystemMatrixSPVectorizedIdentityMPI : public SGPP::datadriven::DMSystemMatrixBaseSP {
      private:
        /// vectorization mode
        VectorizationType vecMode_;
        /// Number of orignal training instances
        size_t numTrainingInstances_;
        /// Number of patched and used training instances
        size_t numPatchedTrainingInstances_;
        /// OperationB for calculating the data matrix
        SGPP::parallel::OperationMultipleEvalVectorizedSP* B_;

        double waitting_time;
      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to SGPP::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         * @param vecMode vectorization mode
         */
        DMSystemMatrixSPVectorizedIdentityMPI(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrixSPVectorizedIdentityMPI();

        virtual void mult(SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& result);

        virtual void generateb(SGPP::base::DataVectorSP& classes, SGPP::base::DataVectorSP& b);

        virtual void rebuildLevelAndIndex();

      private:
        /// how to distribute storage array
        int* _mpi_grid_sizes;
        int* _mpi_grid_offsets;

        /// reference to grid. needed to get new grid size after it changes
        SGPP::base::Grid& m_grid;

        /// how to distribute dataset
        int* _mpi_data_sizes;
        int* _mpi_data_offsets;

        /**
         * Wrapper function that handles communication after calculation and time measurement
         */
        void multVec(base::DataVectorSP& alpha, base::DataVectorSP& result);

        /**
         * Wrapper function that handles communication after calculation and time measurement
         */
        void multTransposeVec(base::DataVectorSP& source, base::DataVectorSP& result);
    };

  }
}

#endif
#endif /* DMSYSTEMMATRIXSPVECTORIZEDIDENTITYMPI_HPP */

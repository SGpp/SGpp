/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITY_HPP
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITY_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/datadriven/algorithm/DMSystemMatrixBaseSP.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <string>

namespace sg {
  namespace parallel {

    /**
     * Class that implements the virtual class sg::base::OperationMatrix for the
     * application of classification for the Systemmatrix
     *
     * The Identity matrix is used as regularization operator.
     *
     * For the Operation B's mult and mutlTransposed functions
     * vectorized formulations are used.
     *
     * In this class single precision DataVectors are used.
     */
    class DMSystemMatrixSPVectorizedIdentity : public sg::datadriven::DMSystemMatrixBaseSP {
      private:
        /// vectorization mode
        VectorizationType vecMode_;
        /// Number of original training instances
        size_t numTrainingInstances_;
        /// Number of patched and used training instances
        size_t numPatchedTrainingInstances_;
        /// OperationB for calculating the data matrix
        sg::parallel::OperationMultipleEvalVectorizedSP* B_;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to sg::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         * @param vecMode vectorization mode, possible values are SSE, AVX, OCL, ArBB
         */
        DMSystemMatrixSPVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrixSPVectorizedIdentity();

        virtual void mult(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result);

        virtual void generateb(sg::base::DataVectorSP& classes, sg::base::DataVectorSP& b);

        virtual void rebuildLevelAndIndex();
    };

  }
}

#endif /* DMSYSTEMMATRIXSPVECTORIZEDIDENTITY_HPP */

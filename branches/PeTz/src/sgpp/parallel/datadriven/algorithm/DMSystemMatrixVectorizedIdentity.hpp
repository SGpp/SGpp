/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP
#define DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/tools/TypesParallel.hpp"

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
     */
    class DMSystemMatrixVectorizedIdentity : public sg::datadriven::DMSystemMatrixBase {
      private:
        /// vectorization mode
        VectorizationType vecMode_;
        /// Number of orignal training instances
        size_t numTrainingInstances_;
        /// Number of patched and used training instances
        size_t numPatchedTrainingInstances_;
        /// OperationB for calculating the data matrix
        sg::parallel::OperationMultipleEvalVectorized* B_;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to sg::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         * @param vecMode vectorization mode
         */
        DMSystemMatrixVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrixVectorizedIdentity();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);

        virtual void rebuildLevelAndIndex();
    };

  }
}

#endif /* DMSYSTEMMATRIXVECTORIZEDIDENTITY_HPP */

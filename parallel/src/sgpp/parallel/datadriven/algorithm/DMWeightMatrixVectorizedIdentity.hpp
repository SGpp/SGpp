/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)

#ifndef DMWEIGHTMATRIXVECTORIZEDIDENTITY_HPP
#define DMWEIGHTMATRIXVECTORIZEDIDENTITY_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <string>

namespace sg {
  namespace parallel {

    /**
     * Class that implements the virtual class sg::base::OperationMatrix for the
     * application of classification for the Systemmatrix with weight
     *
     * The Identity matrix is used as regularization operator.
     *
     * For the Operation B's mult and mutlTransposed functions
     * vectorized formulations in SSE, AVX, OpenCL or Intel Array Building Blocks
     * are used.
     */
    class DMWeightMatrixVectorizedIdentity : public sg::base::OperationMatrix {
      private:
        /// the lambda, the regularisation parameter
        double lamb;
        /// OperationB for calculating the data matrix
        sg::parallel::OperationMultipleEvalVectorized* B;
        /// Pointer to the data vector
        sg::base::DataMatrix* data;
        /// Pointer to the weight vector
        sg::base::DataVector* weight;
        /// Number of orignal training instances
        size_t numTrainingInstances;
        /// Number of patched and used training instances
        size_t numPatchedTrainingInstances;
        /// vectorization mode, possible values are SSE, AVX, OCL, ArBB
        VectorizationType vecMode;
        /// vector width, class internal variable to enable padding and patching of vectors
        size_t vecWidth;
        // save some timings during computation
        /// time needed for Mult
        double completeTimeMult;
        /// time needed only for the computation of mult, interesting on accelerator boards
        double computeTimeMult;
        /// time needed for Mult transposed
        double completeTimeMultTrans;
        /// time needed only for the computation of mult transposed, interesting on accelerator boards
        double computeTimeMultTrans;
        /// Stopwatch needed to determine the durations of mult and mult transposed
        sg::base::SGppStopwatch* myTimer;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to sg::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         * @param w the weights to the training data
         * @param vecMode vectorization mode, possible values are X86SIMD, OCL, ArBB, HYBRID_X86SIMD_OCL
         */
        DMWeightMatrixVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, sg::base::DataVector& w, VectorizationType vecMode);

        /**
         * Std-Destructor
         */
        virtual ~DMWeightMatrixVectorizedIdentity();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);

        /**
         * rebuilds the sg::base::DataMatrix for Level and Index
         */
        void rebuildLevelAndIndex();

        /**
         * resets all timers to 0
         */
        void resetTimers();

        /**
         * gets the timer's values by saving them into call by reference values
         *
         * @param timeMult variable to store overall time needed for Mult
         * @param computeMult variable to store compute time needed for Mult
         * @param timeMultTrans variable to store everall time needed for Mult Transposed
         * @param computeMultTrans variable to store compute time needed for Mult Transposed
         */
        void getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans);
    };

  }
}

#endif /* DMWEIGHTMATRIXVECTORIZEDIDENTITY_HPP */

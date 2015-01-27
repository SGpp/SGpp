// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef DMWEIGHTMATRIXSPVECTORIZEDIDENTITY_HPP
#define DMWEIGHTMATRIXSPVECTORIZEDIDENTITY_HPP

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrixSP.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <string>

#include <sgpp/globaldef.hpp>


// TODO David
#if USE_DOUBLE_PRECISION==0

namespace SGPP {
  namespace parallel {

    /**
     * Class that implements the virtual class SGPP::base::OperationMatrix for the
     * application of classification for the Systemmatrix with weight
     *
     * The Identity matrix is used as regularization operator.
     *
     * For the Operation B's mult and mutlTransposed functions
     * vectorized formulations in SSE, AVX, OpenCL or Intel Array Building Blocks
     * are used.
     *
     * Sinple precision variant!
     */
    class DMWeightMatrixSPVectorizedIdentity : public SGPP::base::OperationMatrixSP {
      private:
        /// the lambda, the regularisation parameter
        float lamb;
        /// OperationB for calculating the data matrix
        SGPP::parallel::OperationMultipleEvalVectorizedSP* B;
        /// Pointer to the data vector
        SGPP::base::DataMatrixSP* data;
        /// Pointer to the weight vector
        SGPP::base::DataVectorSP* weight;
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
        SGPP::base::SGppStopwatch* myTimer;

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param trainData reference to SGPP::base::DataMatrix that contains the training data
         * @param lambda the lambda, the regression parameter
         * @param w the weights to the training data
         * @param vecMode vectorization mode, possible values are X86SIMD, OCL, ArBB, HYBRID_X86SIMD_OCL
         */
        DMWeightMatrixSPVectorizedIdentity(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrixSP& trainData, float lambda, SGPP::base::DataVectorSP& w, VectorizationType vecMode);

        /**
         * Std-Destructor
         */
        virtual ~DMWeightMatrixSPVectorizedIdentity();

        virtual void mult(SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        void generateb(SGPP::base::DataVectorSP& classes, SGPP::base::DataVectorSP& b);

        /**
         * rebuilds the SGPP::base::DataMatrix for Level and Index
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

#endif
#endif /* DMWEIGHTMATRIXSPVECTORIZEDIDENTITY_HPP */

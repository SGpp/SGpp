// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMSYSTEMMATRIXBASE_HPP
#define DMSYSTEMMATRIXBASE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * Abstract class that defines the virtual class SGPP::base::OperationMatrix for
     * classification and regression problems
     */
    class DMSystemMatrixBase : public SGPP::base::OperationMatrix {
      protected:
        /// the dataset
        SGPP::base::DataMatrix* dataset_;
        /// the lambda, the regularisation parameter
        float_t lambda_;
        /// time needed for Mult
        float_t completeTimeMult_;
        /// time needed only for the computation of mult, interesting on accelerator boards
        float_t computeTimeMult_;
        /// time needed for Mult transposed
        float_t completeTimeMultTrans_;
        /// time needed only for the computation of mult transposed, interesting on accelerator boards
        float_t computeTimeMultTrans_;
        /// Stopwatch needed to determine the durations of mult and mult transposed
        SGPP::base::SGppStopwatch* myTimer_;

      public:
        /**
         * Std-Constructor
         *
         * @param trainData matrix with training data
         * @param lambda the lambda, the regression parameter
         */
        DMSystemMatrixBase(SGPP::base::DataMatrix& trainData, float_t lambda);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrixBase();

        virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        virtual void generateb(SGPP::base::DataVector& classes, SGPP::base::DataVector& b) = 0;

        /**
         * forward declaration
         *
         * rebuilds the SGPP::base::DataMatrix for Level and Index
         * this routine is needed for supporting adaptiva grids
         * with vectorized high performance kernels
         */
        virtual void prepareGrid();

        /**
         * resets all timers to 0
         */
        virtual void resetTimers();

        /**
         * gets the timer's values by saving them into call by reference values
         *
         * @param timeMult variable to store overall time needed for Mult
         * @param computeMult variable to store compute time needed for Mult
         * @param timeMultTrans variable to store everall time needed for Mult Transposed
         * @param computeMultTrans variable to store compute time needed for Mult Transposed
         */
        virtual void getTimers(float_t& timeMult, float_t& computeMult, float_t& timeMultTrans, float_t& computeMultTrans);
    };

  }
}

#endif /* DMSYSTEMMATRIXBASE_HPP */

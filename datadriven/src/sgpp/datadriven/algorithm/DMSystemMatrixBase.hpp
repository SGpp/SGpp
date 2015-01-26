/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXBASE_HPP
#define DMSYSTEMMATRIXBASE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

namespace sg {
  namespace datadriven {

    /**
     * Abstract class that defines the virtual class sg::base::OperationMatrix for
     * classification and regression problems
     */
    class DMSystemMatrixBase : public sg::base::OperationMatrix {
      protected:
        /// the dataset
        sg::base::DataMatrix *dataset_;
        /// the lambda, the regularisation parameter
        double lambda_;
        /// time needed for Mult
        double completeTimeMult_;
        /// time needed only for the computation of mult, interesting on accelerator boards
        double computeTimeMult_;
        /// time needed for Mult transposed
        double completeTimeMultTrans_;
        /// time needed only for the computation of mult transposed, interesting on accelerator boards
        double computeTimeMultTrans_;
        /// Stopwatch needed to determine the durations of mult and mult transposed
        sg::base::SGppStopwatch* myTimer_;

      public:
        /**
         * Std-Constructor
         *
         * @param trainData matrix with training data
         * @param lambda the lambda, the regression parameter
         */
        DMSystemMatrixBase(sg::base::DataMatrix& trainData, double lambda);

        /**
         * Std-Destructor
         */
        virtual ~DMSystemMatrixBase();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b) = 0;

        /**
         * forward declaration
         *
         * rebuilds the sg::base::DataMatrix for Level and Index
         * this routine is needed for supporting adaptiva grids
         * with vectorized high performance kernels
         */
        virtual void rebuildLevelAndIndex();

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
        virtual void getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans);
    };

  }
}

#endif /* DMSYSTEMMATRIXBASE_HPP */

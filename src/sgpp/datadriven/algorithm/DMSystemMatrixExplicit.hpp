/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)


#ifndef DMSYSTEMMATRIXEXPLICIT_HPP
#define DMSYSTEMMATRIXEXPLICIT_HPP

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"
#include "base/grid/GridStorage.hpp"
#include <limits.h>

namespace sg {
  namespace datadriven {

    class DMSystemMatrixExplicit: public sg::datadriven::DMSystemMatrixBase {
      public:
        /**
         * Std-Constructor
         *
         * @param lambda the lambda, the regression parameter
         */
        DMSystemMatrixExplicit(sg::base::DataMatrix& trainData, double lambda);

        /**
         * Std-Destructor
         */
        virtual  ~DMSystemMatrixExplicit();

        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

        /**
         * Generates the right hand side of the classification equation
         *
         * @param classes the class information of the training data
         * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
         */
        virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b);

        virtual double* compileMatrix(sg::base::GridStorage* storage,unsigned int first_index = 0, unsigned int last_index = UINT_MAX);

        /*virtual void serialize(std::ostream& ostr);
        virtual void serialize(const std::string& ostr);

        virtual void unserialize(std::ostream& ostr);
        virtual void unserialize(const std::string& ostr);*/

    };

  }
}

#endif /* DMSYSTEMMATRIXEXPLICIT_HPP */

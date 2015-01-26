/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMATRIX_HPP
#define OPERATIONMATRIX_HPP

#include "base/datatypes/DataVector.hpp"

namespace sg {
  namespace base {

    /**
     * Abstract definition of a matrix operator interface.
     * Every time you need to apply a matrix to the ansatzfunction's
     * coefficients derive a class from OperationMatrix
     */
    class OperationMatrix {
      public:
        /**
         * Constructor
         */
        OperationMatrix() {}

        /**
         * Destructor
         */
        virtual ~OperationMatrix() {}

        /**
         * starts the Multiplication with the matrix
         *
         * @param alpha DataVector that contains the ansatzfunctions' coefficients
         * @param result DataVector into which the result of the Laplace operation is stored
         */
        virtual void mult(DataVector& alpha, DataVector& result) = 0;
    };

  }
}

#endif /* OPERATIONMATRIX_HPP */

/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMATRIXSP_HPP
#define OPERATIONMATRIXSP_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>

namespace sg {
  namespace base {

    /**
     * Abstract definition of a matrix operator interface.
     * Everytime you need to apply a matrix to the ansatzfunction's
     * coefficients in single precision derive a class from OperationMatrixSP.
     *
     * This is an re-implementation of the standard OperationMatrix
     * for single precision floating point numbers in order to
     * increase support for GPUs.
     */
    class OperationMatrixSP {
      public:
        /**
         * Constructor
         */
        OperationMatrixSP() {}

        /**
         * Destructor
         */
        virtual ~OperationMatrixSP() {}

        /**
         * starts the Multiplication with the Laplace matrix
         *
         * @param alpha DataVectorSP that contains the ansatzfunctions' coefficients
         * @param result DataVectorSP into which the result of the Laplace operation is stored
         */
        virtual void mult(DataVectorSP& alpha, DataVectorSP& result) = 0;
    };

  }
}

#endif /* OPERATIONMATRIXSP_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMATRIXSP_HPP
#define OPERATIONMATRIXSP_HPP

#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
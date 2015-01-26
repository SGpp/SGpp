/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARYOCL_HPP
#define OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARYOCL_HPP

#include <sgpp/base/operation/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>

namespace sg {
  namespace parallel {

    /**
     * Implements the standard L 2 scalar product on linear boundary grids
     *
     * @version $HEAD$
     */
    class OperationLTwoDotProductVectorizedLinearBoundaryOCL: public sg::base::OperationMatrix {

      private:
        sg::base::GridStorage* storage;
        sg::base::DataMatrix* level_;
        sg::base::DataMatrix* level_int_;
        sg::base::DataMatrix* index_;
        double* lcl_q;
        OCLPDEKernels OCLPDEKernelsHandle;

        void mult_dirichlet(sg::base::DataVector& alpha, sg::base::DataVector& result);

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationLTwoDotProductVectorizedLinearBoundaryOCL(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLTwoDotProductVectorizedLinearBoundaryOCL();

      protected:
        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);
    };

  }
}

#endif /* OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARYOCL_HPP */

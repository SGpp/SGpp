/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Jacob Jepsen (jepsen@diku.dk)

#ifndef OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAROCL_HPP
#define OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAROCL_HPP

#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/pde/basis/common/OCLPDEKernels.hpp"

namespace sg {
  namespace parallel {

    /**
     * Implements the standard L 2 scalar product on linear grids (no boundaries)
     *
     * @version $HEAD$
     */
    class OperationLTwoDotProductVectorizedLinearOCL : public sg::base::OperationMatrix {

      private:
        sg::base::GridStorage* storage;
        sg::base::DataMatrix* level_;
        sg::base::DataMatrix* level_int_;
        sg::base::DataMatrix* index_;
        double* lcl_q;
        OCLPDEKernels OCLPDEKernelsHandle ;

      public:
        /**
         * Constructor
         *
         * @param storage the grid's sg::base::GridStorage object
         */
        OperationLTwoDotProductVectorizedLinearOCL(sg::base::GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~OperationLTwoDotProductVectorizedLinearOCL();

      protected:
        virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);


    };

  }
}

#endif /* OPERATIONLTWODOTPRODUCTVECTORIZEDLINEAROCL_HPP */

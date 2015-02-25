// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONROSENBLATTTRANSFORMATIONLINEAR_HPP
#define OPERATIONROSENBLATTTRANSFORMATIONLINEAR_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/OperationRosenblattTransformation.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /**
     * keep applying marginalize to function until it's reduced to only 1 dimension
     */

    class OperationRosenblattTransformationLinear: public OperationRosenblattTransformation {
      public:
        OperationRosenblattTransformationLinear(base::Grid* grid) :
          grid(grid) {
        }
        virtual ~OperationRosenblattTransformationLinear() {
        }

        /**
         * Transformation with mixed starting dimensions
         *
         * @param alpha Coefficient vector for current grid
           * @param points Input Matrix
           * @param pointscdf Output Matrix
         */
        void doTransformation(base::DataVector* alpha, base::DataMatrix* points,
                              base::DataMatrix* pointscdf);

        /**
         * Transformation with specified starting dimension
         *
         * @param alpha Coefficient vector for current grid
         * @param points Input Matrix
         * @param pointscdf Output Matrix
         * @param dim_start Starting dimension
         */
        void doTransformation(base::DataVector* alpha, base::DataMatrix* points,
                              base::DataMatrix* pointscdf, size_t dim_start);

      protected:
        base::Grid* grid;
        void doTransformation_start_dimX(base::Grid* g_in, base::DataVector* a_in,
                                         size_t dim_start, base::DataVector* coords1d,
                                         base::DataVector* cdfs1d);
        void doTransformation_in_next_dim(base::Grid* g_in, base::DataVector* a_in,
                                          size_t dim_x, base::DataVector* coords1d, base::DataVector* cdfs1d,
                                          size_t& curr_dim);
        virtual float_t doTransformation1D(base::Grid* grid1d,
                                           base::DataVector* alpha1d, float_t coord1d);
    };

  }
}
#endif /* OPERATIONROSENBLATTTRANSFORMATIONLINEAR_HPP */

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONINVERSEROSENBLATTTRANSFORMATIONLINEAR1D_HPP
#define OPERATIONINVERSEROSENBLATTTRANSFORMATIONLINEAR1D_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/OperationTransformation1D.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace datadriven {

    class OperationInverseRosenblattTransformation1DLinear : public OperationTransformation1D {
      protected:
        base::Grid* grid;
      public:
        OperationInverseRosenblattTransformation1DLinear(base::Grid* grid);
        virtual ~OperationInverseRosenblattTransformation1DLinear();

        /**
         * Inverse Rosenblatt Transformation 1D
         * @param alpha1d
         * @param coord1d
         * @return
         */
        float_t doTransformation1D(base::DataVector* alpha1d, float_t coord1d);
    };

  } /* namespace datadriven */
} /* namespace SGPP */

#endif /* OPERATIONINVERSEROSENBLATTTRANSFORMATIONLINEAR1D_HPP */

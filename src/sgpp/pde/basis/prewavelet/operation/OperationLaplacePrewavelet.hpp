/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACEPREWAVELET_HPP
#define OPERATIONLAPLACEPREWAVELET_HPP

#include "pde/algorithm/UpDownOneOpDimWithShadow.hpp"

#include "base/operation/OperationMatrix.hpp"

#include "base/algorithm/sweep.hpp"

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

#include<iostream>

using namespace sg::base;
namespace sg {
  namespace pde {

    /**
     * Implementation for linear functions of Laplace Operation, prewavelet grids without boundaries.
     * With prewavelets the calculation of the gradient part of the up down algorithm is the more complicated
     * one whereas the normal part is eased. For details on the implementation please refer to the documentation
     * of the detail-classes LaplaceDownGradientPrewavelet.hpp, LaplaceUpGradientPrewavelet.hpp and
     * LaplaceDownPrewavelet.hpp.
     */
    class OperationLaplacePrewavelet : public UpDownOneOpDimWithShadow {
      public:
        /**
         * Constructor of OperationLaplacePrewavelet
         *
         * @param storage Pointer to the grid's gridstorage obejct
         * @param shadowstorage shadow storage fuer prewavelets
         */
        OperationLaplacePrewavelet(sg::base::GridStorage* storage, sg::base::GridStorage* shadowstorage);

        /**
         * Destructor
         */
        virtual ~OperationLaplacePrewavelet();

      protected:

        virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

        virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
    };

  }
}

#endif /* OPERATIONLAPLACEPREWAVELET_HPP */
